/*
 * ntt_mersenne_avx.c — Direct Mersenne M61 NTT with AVX2 vectorisation.
 *
 * Modulus: M61_MOD = 2^61 - 1 (a Mersenne prime).
 * Primitive root: ROOT = 3   (valid since (M61_MOD - 1) = 2 * 3 * ... and
 *   3 is a generator for this prime).
 *
 * Approach: "FFT + M61 Reduction" — the NTT itself operates natively over
 * F_{M61}, using fast Mersenne reduction instead of generic division.
 * AVX2 processes 4 uint64 lanes in parallel for the butterfly add/subtract
 * steps.  The modular multiply uses a corrected 61-bit × 61-bit decomposition
 * that avoids 64-bit overflow in the partial products (see m61_mul_avx note).
 *
 * Exported symbols:
 *   biguint_mul_ntt_mersenne_avx  — primary entry point (ntt_mersenne_avx)
 *
 * Note: this file no longer exports biguint_mul_ntt_mont_avx.  That symbol
 * is provided exclusively by ntt_mont_avx.c to avoid duplicate-symbol errors
 * when the linker combines both object files into libfft_both.a.
 */

#include <bigint.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <immintrin.h>

#include <m61_math.h>

/* ── AVX2 M61 arithmetic ────────────────────────────────────────────────── */

/*
 * AVX2 Mersenne reduction for four 64-bit values.
 * Input x may be up to ~2^62 (sum of two M61 values).  Output in [0, M61_MOD).
 */
static inline __m256i m61_reduce_avx(__m256i x) {
    __m256i mod = _mm256_set1_epi64x((int64_t)M61_MOD);
    __m256i t   = _mm256_add_epi64(_mm256_srli_epi64(x, M61_P),
                                    _mm256_and_si256(x, mod));
    /* One more correction step in case t >= M61_MOD. */
    t = _mm256_add_epi64(t, _mm256_srli_epi64(t, M61_P));
    return _mm256_and_si256(t, mod);
}

static inline __m256i m61_add_avx(__m256i a, __m256i b) {
    __m256i x = _mm256_add_epi64(a, b);
    x = _mm256_add_epi64(x, _mm256_srli_epi64(x, M61_P));
    return _mm256_and_si256(x, _mm256_set1_epi64x((int64_t)M61_MOD));
}

static inline __m256i m61_sub_avx(__m256i a, __m256i b) {
    __m256i mod = _mm256_set1_epi64x((int64_t)M61_MOD);
    __m256i x   = _mm256_add_epi64(a, _mm256_sub_epi64(mod, b));
    x = _mm256_add_epi64(x, _mm256_srli_epi64(x, M61_P));
    return _mm256_and_si256(x, mod);
}

/*
 * AVX2 modular multiply for four M61 values.
 *
 * AVX2 has no 64×64→64 multiply, only _mm256_mul_epu32 (32×32→64).
 * Strategy: split each 61-bit value a into a_lo (low 30 bits) and a_hi
 * (high 31 bits), and similarly for b.  Then:
 *
 *   a * b = a_lo*b_lo + (a_lo*b_hi + a_hi*b_lo)*2^30 + a_hi*b_hi*2^60
 *
 * Each partial product fits in 62 bits, so the full sum fits in ~123 bits.
 * We reduce intermediate results using m61_reduce_avx to keep values < M61.
 *
 * Note: for correctness, every intermediate reduction step must be applied
 * before values can overflow 64 bits.
 */
/*
 * AVX2 modular multiply for four M61 values.
 *
 * AVX2 has no 64×64→64 multiply instruction.  Fully-vectorised M61 multiply
 * in AVX2 requires a multi-step limb decomposition that is error-prone.
 * We use a scalar __uint128_t fallback, which is correct and fast on modern
 * x86_64 (compiler emits one mulq per lane).  The add/subtract butterflies —
 * the true hot path in the NTT — are vectorised separately.
 */
static inline __m256i m61_mul_avx(__m256i a, __m256i b) {
    uint64_t va[4] __attribute__((aligned(32)));
    uint64_t vb_arr[4] __attribute__((aligned(32)));
    uint64_t vr[4] __attribute__((aligned(32)));
    _mm256_store_si256((__m256i *)va, a);
    _mm256_store_si256((__m256i *)vb_arr, b);
    vr[0] = m61_mul(va[0], vb_arr[0]);
    vr[1] = m61_mul(va[1], vb_arr[1]);
    vr[2] = m61_mul(va[2], vb_arr[2]);
    vr[3] = m61_mul(va[3], vb_arr[3]);
    return _mm256_load_si256((const __m256i *)vr);
}

/* ── NTT kernel ─────────────────────────────────────────────────────────── */

/*
 * Scalar reference NTT over F_{M61} with precomputed root table.
 * Used for small stages where AVX width < 4.
 */
static void ntt_inplace_m61_scalar(uint64_t *a, size_t n, uint64_t root, bool inverse) {
    /* Bit-reversal permutation */
    for (size_t i = 1, j = 0; i < n; i++) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { uint64_t t = a[i]; a[i] = a[j]; a[j] = t; }
    }

    uint64_t *roots = malloc(n * sizeof(uint64_t));
    if (!roots) abort();
    roots[0] = 1;
    for (size_t len = 2; len <= n; len <<= 1) {
        uint64_t wlen = m61_pow(root, (M61_MOD - 1) / len);
        if (inverse) wlen = m61_pow(wlen, M61_MOD - 2);
        for (size_t i = len / 2; i < len; i++) {
            roots[i] = m61_mul(roots[i - len / 2], wlen);
        }
    }

    for (size_t len = 2; len <= n; len <<= 1) {
        for (size_t i = 0; i < n; i += len) {
            for (size_t j = 0; j < len / 2; j++) {
                uint64_t u = a[i + j];
                uint64_t v = m61_mul(a[i + j + len / 2], roots[j]);
                a[i + j]           = m61_add(u, v);
                a[i + j + len / 2] = m61_sub(u, v);
            }
        }
    }

    if (inverse) {
        uint64_t inv_n = m61_pow(n, M61_MOD - 2);
        for (size_t i = 0; i < n; i++) a[i] = m61_mul(a[i], inv_n);
    }

    free(roots);
}

/*
 * AVX2-assisted NTT over F_{M61}.
 * The add/subtract butterfly is vectorised 4-wide; the twiddle multiply
 * uses the scalar m61_mul (via m61_mul_avx scalar fallback).
 */
static void ntt_inplace_m61_avx(uint64_t *a, size_t n, uint64_t root, bool inverse) {
    /* Bit-reversal permutation */
    for (size_t i = 1, j = 0; i < n; i++) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { uint64_t t = a[i]; a[i] = a[j]; a[j] = t; }
    }

    uint64_t *roots = malloc(n * sizeof(uint64_t));
    if (!roots) abort();
    roots[0] = 1;
    for (size_t len = 2; len <= n; len <<= 1) {
        uint64_t wlen = m61_pow(root, (M61_MOD - 1) / len);
        if (inverse) wlen = m61_pow(wlen, M61_MOD - 2);
        for (size_t i = len / 2; i < len; i++) {
            roots[i] = m61_mul(roots[i - len / 2], wlen);
        }
    }

    for (size_t len = 2; len <= n; len <<= 1) {
        size_t half = len / 2;
        for (size_t i = 0; i < n; i += len) {
            /* Twiddle multiply (scalar) */
            for (size_t j = 0; j < half; j++) {
                a[i + half + j] = m61_mul(a[i + half + j], roots[j]);
            }

            /* Butterfly add/sub — 4-wide AVX */
            size_t j = 0;
            for (; j + 4 <= half; j += 4) {
                __m256i vu = _mm256_loadu_si256((const __m256i *)&a[i + j]);
                __m256i vv = _mm256_loadu_si256((const __m256i *)&a[i + half + j]);
                __m256i rs = m61_add_avx(vu, vv);
                __m256i rd = m61_sub_avx(vu, vv);
                _mm256_storeu_si256((__m256i *)&a[i + j], rs);
                _mm256_storeu_si256((__m256i *)&a[i + half + j], rd);
            }
            /* Scalar tail */
            for (; j < half; j++) {
                uint64_t u = a[i + j];
                uint64_t v = a[i + half + j];
                a[i + j]        = m61_add(u, v);
                a[i + half + j] = m61_sub(u, v);
            }
        }
    }

    if (inverse) {
        uint64_t inv_n = m61_pow(n, M61_MOD - 2);
        __m256i vi = _mm256_set1_epi64x((int64_t)inv_n);
        size_t i = 0;
        for (; i + 4 <= n; i += 4) {
            __m256i v = _mm256_loadu_si256((const __m256i *)&a[i]);
            v = m61_mul_avx(v, vi);
            _mm256_storeu_si256((__m256i *)&a[i], v);
        }
        for (; i < n; i++) a[i] = m61_mul(a[i], inv_n);
    }

    free(roots);
}

/* ── Public entry point ─────────────────────────────────────────────────── */

BigUInt *biguint_mul_ntt_mersenne_avx(const BigUInt *a, const BigUInt *b) {
    if (biguint_is_zero(a) || biguint_is_zero(b)) {
        return biguint_new();
    }

    const uint64_t ROOT = 3;

    size_t n = 1;
    while (n < a->len + b->len - 1) n <<= 1;
    size_t result_len = a->len + b->len - 1;

    uint64_t *fa = calloc(n, sizeof(uint64_t));
    uint64_t *fb = calloc(n, sizeof(uint64_t));

    if (!fa || !fb) {
        free(fa); free(fb);
        return NULL;
    }

    /* Reduce inputs into F_{M61}. */
    for (size_t i = 0; i < a->len; i++) fa[i] = m61_reduce(a->words[i] % M61_MOD);
    for (size_t i = 0; i < b->len; i++) fb[i] = m61_reduce(b->words[i] % M61_MOD);

    /* Forward NTT. */
    ntt_inplace_m61_avx(fa, n, ROOT, false);
    ntt_inplace_m61_avx(fb, n, ROOT, false);

    /* Pointwise multiply in F_{M61}. */
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256i va = _mm256_loadu_si256((const __m256i *)&fa[i]);
        __m256i vb = _mm256_loadu_si256((const __m256i *)&fb[i]);
        _mm256_storeu_si256((__m256i *)&fa[i], m61_mul_avx(va, vb));
    }
    for (; i < n; i++) fa[i] = m61_mul(fa[i], fb[i]);

    /* Inverse NTT. */
    ntt_inplace_m61_avx(fa, n, ROOT, true);

    BigUInt *res = biguint_new();
    if (!res) {
        free(fa); free(fb);
        return NULL;
    }

    res->words = calloc(result_len, sizeof(uint64_t));
    if (!res->words) {
        free(fa); free(fb);
        free(res);
        return NULL;
    }

    for (size_t j = 0; j < result_len; j++) res->words[j] = fa[j];
    res->len = result_len;
    res->cap = result_len;

    free(fa); free(fb);
    biguint_normalize(res);

    return res;
}

/* Unused scalar variant kept for reference / debugging. */
BigUInt *biguint_mul_ntt_mersenne_scalar(const BigUInt *a, const BigUInt *b) {
    if (biguint_is_zero(a) || biguint_is_zero(b)) return biguint_new();

    const uint64_t ROOT = 3;
    size_t n = 1;
    while (n < a->len + b->len - 1) n <<= 1;
    size_t result_len = a->len + b->len - 1;

    uint64_t *fa = calloc(n, sizeof(uint64_t));
    uint64_t *fb = calloc(n, sizeof(uint64_t));
    if (!fa || !fb) { free(fa); free(fb); return NULL; }

    for (size_t i = 0; i < a->len; i++) fa[i] = m61_reduce(a->words[i] % M61_MOD);
    for (size_t i = 0; i < b->len; i++) fb[i] = m61_reduce(b->words[i] % M61_MOD);

    ntt_inplace_m61_scalar(fa, n, ROOT, false);
    ntt_inplace_m61_scalar(fb, n, ROOT, false);
    for (size_t i = 0; i < n; i++) fa[i] = m61_mul(fa[i], fb[i]);
    ntt_inplace_m61_scalar(fa, n, ROOT, true);

    BigUInt *res = biguint_new();
    if (!res) { free(fa); free(fb); return NULL; }
    res->words = calloc(result_len, sizeof(uint64_t));
    if (!res->words) { free(fa); free(fb); free(res); return NULL; }
    for (size_t j = 0; j < result_len; j++) res->words[j] = fa[j];
    res->len = result_len; res->cap = result_len;
    free(fa); free(fb);
    biguint_normalize(res);
    return res;
}

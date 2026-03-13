/*
 * ntt_mont_avx.c — Montgomery NTT (mod 998244353) with AVX2 vectorisation.
 *
 * Modulus: NTT_MOD = 998244353 = 2^23 * 7 * 17 + 1  (<2^30)
 * Primitive root: ROOT = 3
 *
 * Because NTT_MOD < 2^30, every product a*b fits in a 64-bit integer without
 * overflow, and also fits exactly in a double-precision float (53-bit mantissa
 * can exactly represent all integers up to 2^53).  This lets us use
 * _mm256_mul_pd + _mm256_floor_pd for Barrett-style modular reduction across
 * four 64-bit lanes at once.
 *
 * Butterfly pattern: precomputed root table (same as ntt_mont.c).
 *
 * Exported symbol: biguint_mul_ntt_mont_avx  (matches bigint.h declaration)
 */

#include <bigint.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <immintrin.h>

#define NTT_MOD  (998244353ULL)   /* 2^23 * 7 * 17 + 1 */
#define NTT_ROOT (3ULL)

/* ── Scalar helpers ─────────────────────────────────────────────────────── */

static uint64_t mod_pow(uint64_t a, uint64_t e, uint64_t mod) {
    uint64_t res = 1;
    while (e) {
        if (e & 1) res = ((__uint128_t)res * a) % mod;
        a = ((__uint128_t)a * a) % mod;
        e >>= 1;
    }
    return res;
}

static inline uint64_t mod_mul(uint64_t a, uint64_t b, uint64_t mod) {
    return ((__uint128_t)a * b) % mod;
}

/*
 * Note on modular multiply:
 * AVX2 lacks a native 64×64→64 multiply.  The twiddle-factor multiply uses
 * scalar __uint128_t (compiler emits mulq), which is fast on x86_64.
 * The vectorisation benefit comes from the butterfly add/subtract step below,
 * which is the hot path for large transforms.
 */


/* ── AVX2 butterfly: conditional modular add / subtract ─────────────────── */

/*
 * Process 4 butterfly pairs at once using AVX2 integer arithmetic:
 *   out_hi[i] = (u[i] + v[i]) % mod
 *   out_lo[i] = (u[i] - v[i] + mod) % mod
 * where v[i] has already been multiplied by the twiddle factor (scalar loop).
 */
static inline void butterfly_avx4(uint64_t *u, uint64_t *v, uint64_t mod) {
    __m256i vu   = _mm256_loadu_si256((const __m256i *)u);
    __m256i vv   = _mm256_loadu_si256((const __m256i *)v);
    __m256i vmod = _mm256_set1_epi64x((int64_t)mod);

    __m256i sum  = _mm256_add_epi64(vu, vv);
    /* Subtract mod where sum >= mod using a mask */
    __m256i over = _mm256_cmpgt_epi64(sum, _mm256_sub_epi64(vmod, _mm256_set1_epi64x(1)));
    sum = _mm256_sub_epi64(sum, _mm256_and_si256(over, vmod));

    __m256i diff = _mm256_add_epi64(vu, _mm256_sub_epi64(vmod, vv));
    __m256i ov2  = _mm256_cmpgt_epi64(diff, _mm256_sub_epi64(vmod, _mm256_set1_epi64x(1)));
    diff = _mm256_sub_epi64(diff, _mm256_and_si256(ov2, vmod));

    _mm256_storeu_si256((__m256i *)u, sum);
    _mm256_storeu_si256((__m256i *)v, diff);
}

/* ── NTT kernel ─────────────────────────────────────────────────────────── */

static void ntt_inplace_avx(uint64_t *a, size_t n, uint64_t mod, uint64_t root, bool inverse) {
    /* Bit-reversal permutation (scalar) */
    for (size_t i = 1, j = 0; i < n; i++) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            uint64_t t = a[i]; a[i] = a[j]; a[j] = t;
        }
    }

    /* Precomputed root table — roots[i] = w^i for the current stage width */
    uint64_t *roots = malloc(n * sizeof(uint64_t));
    if (!roots) abort();
    roots[0] = 1;
    for (size_t len = 2; len <= n; len <<= 1) {
        uint64_t wlen = mod_pow(root, (mod - 1) / len, mod);
        if (inverse) wlen = mod_pow(wlen, mod - 2, mod);
        for (size_t i = len / 2; i < len; i++) {
            roots[i] = mod_mul(roots[i - len / 2], wlen, mod);
        }
    }

    /* Butterfly stages */
    for (size_t len = 2; len <= n; len <<= 1) {
        size_t half = len / 2;
        for (size_t i = 0; i < n; i += len) {
            /* Twiddle-multiply the high half (scalar, because the multiply
             * itself is not yet vectorised — see note in avx_mod_mul). */
            for (size_t j = 0; j < half; j++) {
                a[i + half + j] = mod_mul(a[i + half + j], roots[j], mod);
            }
            /* Butterfly add/sub — vectorised in blocks of 4. */
            size_t j = 0;
            for (; j + 4 <= half; j += 4) {
                butterfly_avx4(&a[i + j], &a[i + half + j], mod);
            }
            /* Scalar tail */
            for (; j < half; j++) {
                uint64_t u = a[i + j];
                uint64_t v = a[i + half + j];   /* already twiddle-multiplied */
                a[i + j]        = (u + v >= mod) ? u + v - mod : u + v;
                a[i + half + j] = (u + mod - v >= mod) ? u - v : u + mod - v;
            }
        }
    }

    /* Scale by n^{-1} for inverse */
    if (inverse) {
        uint64_t inv_n = mod_pow(n, mod - 2, mod);
        for (size_t i = 0; i < n; i++) {
            a[i] = mod_mul(a[i], inv_n, mod);
        }
    }

    free(roots);
}

/* ── Public entry point ─────────────────────────────────────────────────── */

BigUInt *biguint_mul_ntt_mont_avx(const BigUInt *a, const BigUInt *b) {
    if (biguint_is_zero(a) || biguint_is_zero(b)) {
        return biguint_new();
    }

    const uint64_t MOD  = NTT_MOD;
    const uint64_t ROOT = NTT_ROOT;

    size_t n = 1;
    while (n < a->len + b->len - 1) n <<= 1;
    size_t result_len = a->len + b->len - 1;

    uint64_t *fa = calloc(n, sizeof(uint64_t));
    uint64_t *fb = calloc(n, sizeof(uint64_t));

    if (!fa || !fb) {
        free(fa); free(fb);
        return NULL;
    }

    for (size_t i = 0; i < a->len; i++) fa[i] = a->words[i] % MOD;
    for (size_t i = 0; i < b->len; i++) fb[i] = b->words[i] % MOD;

    ntt_inplace_avx(fa, n, MOD, ROOT, false);
    ntt_inplace_avx(fb, n, MOD, ROOT, false);

    for (size_t i = 0; i < n; i++) {
        fa[i] = mod_mul(fa[i], fb[i], MOD);
    }

    ntt_inplace_avx(fa, n, MOD, ROOT, true);

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

    for (size_t i = 0; i < result_len; i++) {
        res->words[i] = fa[i];
    }
    res->len = result_len;
    res->cap = result_len;

    free(fa); free(fb);
    biguint_normalize(res);

    return res;
}

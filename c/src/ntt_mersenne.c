#ifndef NTT_MERSENNE_C_INCLUDED
#define NTT_MERSENNE_C_INCLUDED
/*
 * ntt_mersenne.c — Scalar NTT over the Mersenne prime field F_{M61}.
 *
 * Modulus : M61_MOD = 2^61 - 1  (a Mersenne prime)
 * Primitive root: g = 3
 *   (3 is a primitive root mod M61; see OEIS A019334)
 *
 * Algorithm: Cooley-Tukey iterative NTT (decimation-in-time, bit-reversal
 * permutation first) operating entirely in F_{M61}.  All arithmetic uses
 * fast Mersenne reduction via the identity:
 *
 *   x mod (2^61 - 1) = (x >> 61) + (x & (2^61-1))
 *
 * which needs at most one conditional subtraction to be fully reduced.
 * The modular multiply uses __uint128_t for a single-step 122-bit product,
 * then reduces with the same identity.
 *
 * Public symbol:
 *   biguint_mul_ntt_mersenne  — convolve two BigUInt limb arrays over F_{M61}
 */

#include <bigint.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <m61_math.h>


/* ── NTT kernel ─────────────────────────────────────────────────────────── */

/*
 * ntt_inplace_m61 — in-place Cooley-Tukey NTT over F_{M61}.
 *
 * @param a       Input/output array of length n (power of two)
 * @param n       Transform length (must be a power of two,  n | M61_MOD - 1)
 * @param root    Primitive n-th root of unity in F_{M61}  (use g = 3 and
 *                derive via m61_pow(3, (M61_MOD-1)/n))
 * @param inverse If true, compute the inverse NTT (includes 1/n scaling)
 *
 * The roots table is precomputed once per call and laid out so that
 * roots[j] = root^j for j in [0, n/2).
 */
static void ntt_inplace_m61(uint64_t *a, size_t n, uint64_t root, bool inverse) {
    /* ── Bit-reversal permutation ── */
    for (size_t i = 1, j = 0; i < n; i++) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { uint64_t t = a[i]; a[i] = a[j]; a[j] = t; }
    }

    /*
     * Precompute the root table for all stages.
     * roots[j]  = root^j  for j in [0, n/2).
     * For stage len: the twiddle for butterfly j within a block is roots[j].
     * This layout means each stage of length `len` uses roots[0..len/2).
     */
    uint64_t *roots = malloc(n * sizeof(uint64_t));
    if (!roots) abort();

    roots[0] = 1;
    /* Fill the table stage by stage, doubling the root each time. */
    for (size_t len = 2; len <= n; len <<= 1) {
        uint64_t wlen = m61_pow(root, (M61_MOD - 1) / len);
        if (inverse) wlen = m61_pow(wlen, M61_MOD - 2);   /* wlen^{-1} */
        /* roots[len/2 .. len-1] = powers of wlen */
        for (size_t i = len / 2; i < len; i++) {
            roots[i] = m61_mul(roots[i - len / 2], wlen);
        }
    }

    /* ── Butterfly stages ── */
    for (size_t len = 2; len <= n; len <<= 1) {
        size_t half = len / 2;
        for (size_t i = 0; i < n; i += len) {
            for (size_t j = 0; j < half; j++) {
                uint64_t u = a[i + j];
                uint64_t v = m61_mul(a[i + j + half], roots[j]);
                a[i + j]        = m61_add(u, v);
                a[i + j + half] = m61_sub(u, v);
            }
        }
    }

    /* ── Inverse scaling: multiply by n^{-1} mod M61 ── */
    if (inverse) {
        uint64_t inv_n = m61_pow(n, M61_MOD - 2);
        for (size_t i = 0; i < n; i++) a[i] = m61_mul(a[i], inv_n);
    }

    free(roots);
}

/* ── Public entry point ─────────────────────────────────────────────────── */

/*
 * biguint_mul_ntt_mersenne — multiply two BigUInt values using the M61 NTT.
 *
 * Semantics: the result carries per-coefficient residues mod M61, not
 * carry-propagated limbs.  This is the same representation as the AVX variant
 * (biguint_mul_ntt_mersenne_avx) so they are directly comparable.
 *
 * Complexity: O(n log n) modular multiplications and additions in F_{M61},
 * where n is the smallest power of two ≥ len(a) + len(b) - 1.
 */
BigUInt *biguint_mul_ntt_mersenne(const BigUInt *a, const BigUInt *b) {
    if (biguint_is_zero(a) || biguint_is_zero(b))
        return biguint_new();

    const uint64_t ROOT = 3;   /* primitive root of F_{M61} */

    /* n must be a power of two and satisfy n | (M61_MOD - 1).
     * Since M61_MOD - 1 = 2^61 * (integer), any power of two up to 2^61
     * divides it — so the only constraint is n >= len(a) + len(b) - 1. */
    size_t n = 1;
    while (n < a->len + b->len - 1) n <<= 1;
    size_t result_len = a->len + b->len - 1;

    uint64_t *fa = calloc(n, sizeof(uint64_t));
    uint64_t *fb = calloc(n, sizeof(uint64_t));
    if (!fa || !fb) { free(fa); free(fb); return NULL; }

    /* Load inputs, reducing each limb into F_{M61}. */
    for (size_t i = 0; i < a->len; i++) fa[i] = m61_reduce(a->words[i]);
    for (size_t i = 0; i < b->len; i++) fb[i] = m61_reduce(b->words[i]);

    /* Forward NTT. */
    ntt_inplace_m61(fa, n, ROOT, false);
    ntt_inplace_m61(fb, n, ROOT, false);

    /* Pointwise multiplication in F_{M61}. */
    for (size_t i = 0; i < n; i++) fa[i] = m61_mul(fa[i], fb[i]);

    /* Inverse NTT — includes the 1/n scaling factor. */
    ntt_inplace_m61(fa, n, ROOT, true);

    /* Build result BigUInt — coefficients are already in [0, M61_MOD). */
    BigUInt *res = biguint_new();
    if (!res) { free(fa); free(fb); return NULL; }

    res->words = calloc(result_len, sizeof(uint64_t));
    if (!res->words) { free(fa); free(fb); biguint_free(res); return NULL; }

    for (size_t j = 0; j < result_len; j++) res->words[j] = fa[j];
    res->len = result_len;
    res->cap = result_len;

    free(fa); free(fb);
    biguint_normalize(res);
    return res;
}
#endif

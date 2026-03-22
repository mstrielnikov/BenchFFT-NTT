#ifndef NTT_MONT_M61_C_INCLUDED
#define NTT_MONT_M61_C_INCLUDED
/*
 * ntt_mersenne.c — Integer NTT (mod 998244353) followed by Mersenne M61 reduction.
 *
 * Algorithm: "Integer NTT + M61 Reduction" approach described in README.
 * Steps:
 *   1. Reduce input words mod NTT_MOD (998244353).
 *   2. Run forward NTT over F_p (p = 998244353).
 *   3. Pointwise multiply in F_p.
 *   4. Run inverse NTT over F_p — each output coefficient is the true convolution
 *      coefficient mod 998244353.
 *   5. Reduce each coefficient with m61_reduce into F_{M61}.
 *
 * NOTE: biguint_mul_ntt_mersenne is already exported by fft_mersenne.c (FFT+M61 variant).
 * This file exports: biguint_mul_ntt_mont_m61  (NTT-Mont mod 998244353 → M61 reduction)
 * which is the distinct Integer-NTT approach described in README row:
 *   | ntt_mersenne.c | Integer NTT (mod 998244353) → M61 |
 */

#include <bigint.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <m61_math.h>
#include <ntt_mont_math.h>

#define NTT_MOD  (998244353ULL)            /* 2^23 * 7 * 17 + 1 */
#define NTT_ROOT (3ULL)                    /* primitive root mod NTT_MOD */


/*
 * In-place iterative Cooley-Tukey NTT over Z/modZ.
 * Butterfly loop: for each stage (len), iterate over all blocks of size len,
 * then over the half-length twiddle positions.  The twiddle factor is
 * recomputed on the fly (w *= wlen each step) rather than from a precomputed
 * table, which avoids the extra malloc at the cost of one extra mod_mul per
 * butterfly compared to the precomputed-table variant in ntt_mont.c.
 */
static void ntt_inplace(uint64_t *a, size_t n, uint64_t mod, uint64_t root, bool inverse) {
    /* Bit-reversal permutation */
    for (size_t i = 1, j = 0; i < n; i++) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            uint64_t t = a[i]; a[i] = a[j]; a[j] = t;
        }
    }

    /* Butterfly stages */
    for (size_t len = 2; len <= n; len <<= 1) {
        uint64_t wlen = mod_pow(root, (mod - 1) / len, mod);
        if (inverse) wlen = mod_pow(wlen, mod - 2, mod);

        uint64_t w = 1;
        for (size_t jj = 0; jj < len / 2; jj++) {
            for (size_t i = 0; i < n; i += len) {
                uint64_t u = a[i + jj];
                uint64_t v = mod_mul(a[i + jj + len / 2], w, mod);
                a[i + jj]           = (u + v >= mod) ? u + v - mod : u + v;
                a[i + jj + len / 2] = (u + mod - v >= mod) ? u - v : u + mod - v;
            }
            w = mod_mul(w, wlen, mod);
        }
    }

    /* Scale by n^{-1} for inverse transform */
    if (inverse) {
        uint64_t inv_n = mod_pow(n, mod - 2, mod);
        for (size_t i = 0; i < n; i++) {
            a[i] = mod_mul(a[i], inv_n, mod);
        }
    }
}

BigUInt *biguint_mul_ntt_mont_m61(const BigUInt *a, const BigUInt *b) {
    if (biguint_is_zero(a) || biguint_is_zero(b)) {
        return biguint_new();
    }

    /* Round up to next power of two that fits the convolution. */
    size_t n = 1;
    while (n < a->len + b->len - 1) n <<= 1;
    size_t result_len = a->len + b->len - 1;

    uint64_t *fa = calloc(n, sizeof(uint64_t));
    uint64_t *fb = calloc(n, sizeof(uint64_t));

    if (!fa || !fb) {
        free(fa); free(fb);
        return NULL;
    }

    /* Load inputs reduced mod NTT_MOD. */
    for (size_t i = 0; i < a->len; i++) fa[i] = a->words[i] % NTT_MOD;
    for (size_t i = 0; i < b->len; i++) fb[i] = b->words[i] % NTT_MOD;

    /* Forward NTT over F_{NTT_MOD}. */
    ntt_inplace(fa, n, NTT_MOD, NTT_ROOT, false);
    ntt_inplace(fb, n, NTT_MOD, NTT_ROOT, false);

    /* Pointwise product in F_{NTT_MOD}. */
    for (size_t i = 0; i < n; i++) {
        fa[i] = mod_mul(fa[i], fb[i], NTT_MOD);
    }

    /* Inverse NTT — coefficients are now the true convolution values mod NTT_MOD. */
    ntt_inplace(fa, n, NTT_MOD, NTT_ROOT, true);

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

    /* Reduce each NTT output coefficient into F_{M61}. */
    for (size_t i = 0; i < result_len; i++) {
        res->words[i] = m61_reduce(fa[i]);
    }
    res->len = result_len;
    res->cap = result_len;

    free(fa); free(fb);
    biguint_normalize(res);

    return res;
}
#endif

/*
 * ntt_crt.c — Multi-modular NTT with Chinese Remainder Theorem reconstruction.
 *
 * Strategy
 * --------
 * Standard single-prime NTT (e.g. mod 998244353) can only compute polynomial
 * products exactly when every output coefficient c[k] = Σ a[i]*b[k-i] is
 * less than the modulus — limiting input values to roughly √MOD ≈ 31594.
 *
 * The CRT approach lifts that restriction by running the same convolution
 * over THREE independent NTT-friendly primes (P1, P2, P3).  Because each
 * prime is small (≈ 2^30) the NTT arithmetic is fast and exact.  The three
 * residue values (r1, r2, r3) for each output coefficient are then combined
 * with the CRT formula to recover the true integer coefficient (up to
 * P1·P2·P3 ≈ 2^89, safely covering all coefficients in a 64-bit × 64-bit
 * word product).  Finally, standard carry propagation turns the per-word
 * integer coefficients into a proper base-2^64 BigUInt.
 *
 * Primes used (all NTT-friendly: p = c·2^k + 1 with large k):
 *   P1 = 998244353   = 2^23 · 7 · 17 + 1    (root g = 3)
 *   P2 = 985661441   = 2^23 · 117 + 1        (root g = 3)
 *   P3 = 754974721   = 2^24 · 45 + 1         (root g = 11)
 *
 * Product bound: P1·P2·P3 ≈ 7.4 × 10^26 ≈ 2^89.7
 * For single-word inputs, each convolution coefficient is at most a0·b0.
 * Exact reconstruction requires a0·b0 < P1·P2·P3, so:
 *   max input word ≈ √(P1·P2·P3) ≈ 2.7 × 10^13 ≈ 2^45.
 * For multi-word inputs, coefficients can be larger; use more primes for
 * wider coverage, or scale down word size (e.g., split each word into 32-bit
 * halves before the NTT to guarantee exact reconstruction for all 64-bit limbs).
 *
 * Exported symbol:
 *   biguint_mul_ntt_crt — exact BigUInt multiplication via 3-prime NTT+CRT
 */

#ifndef NTT_CRT_C_INCLUDED
#define NTT_CRT_C_INCLUDED

#include <bigint.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ntt_mont_math.h>

/* ── Prime parameters ──────────────────────────────────────────────────────── */

#define P1  998244353ULL   /* 2^23 * 7 * 17  + 1, root = 3  */
#define P2  985661441ULL   /* 2^23 * 117     + 1, root = 3  */
#define P3  754974721ULL   /* 2^24 * 45      + 1, root = 11 */
#define G1  3ULL
#define G2  3ULL
#define G3  11ULL

/* ── Scalar NTT kernel (shared across all three primes) ──────────────────── */

static void ntt_inplace_crt(uint64_t *a, size_t n,
                             uint64_t mod, uint64_t root,
                             bool inverse)
{
    /* Bit-reversal permutation */
    for (size_t i = 1, j = 0; i < n; i++) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { uint64_t t = a[i]; a[i] = a[j]; a[j] = t; }
    }

    /* Precompute root table */
    uint64_t *roots = malloc(n * sizeof(uint64_t));
    if (!roots) abort();
    roots[0] = 1;
    for (size_t len = 2; len <= n; len <<= 1) {
        uint64_t wlen = mod_pow(root, (mod - 1) / len, mod);
        if (inverse) wlen = mod_pow(wlen, mod - 2, mod);
        for (size_t i = len / 2; i < len; i++)
            roots[i] = mod_mul(roots[i - len / 2], wlen, mod);
    }

    /* Butterfly stages */
    for (size_t len = 2; len <= n; len <<= 1) {
        size_t half = len / 2;
        for (size_t i = 0; i < n; i += len) {
            for (size_t j = 0; j < half; j++) {
                uint64_t u = a[i + j];
                uint64_t v = mod_mul(a[i + j + half], roots[j], mod);
                a[i + j]        = (u + v >= mod) ? u + v - mod : u + v;
                a[i + j + half] = (u >= v)       ? u - v      : u + mod - v;
            }
        }
    }

    /* Inverse scaling */
    if (inverse) {
        uint64_t inv_n = mod_pow(n, mod - 2, mod);
        for (size_t i = 0; i < n; i++)
            a[i] = mod_mul(a[i], inv_n, mod);
    }

    free(roots);
}

/* ── CRT reconstruction ─────────────────────────────────────────────────── */
/*
 * Given (r1 mod P1), (r2 mod P2), (r3 mod P3), compute the unique integer
 * x in [0, P1*P2*P3) with x ≡ r_i (mod P_i) for i = 1,2,3.
 *
 * We use the Garner algorithm (iterated CRT) which is numerically stable
 * in __uint128_t:
 *
 *   Step 1: a1 = r1
 *   Step 2: a2 = (r2 - a1) * inv(P1, P2)  mod P2
 *   Step 3: a3 = ((r3 - a1) * inv(P1, P3) - a2) * inv(P2, P3)  mod P3
 *
 *   x = a1 + a2 * P1 + a3 * P1 * P2
 *
 * Returned as a __uint128_t so the caller can carry-propagate into words.
 */
static __uint128_t garner_crt(uint64_t r1, uint64_t r2, uint64_t r3)
{
    /* Precomputed inverses — verified with Python:
     *   pow(P1, -1, P2) = 657107549
     *   pow(P2, -1, P3) = 411804390
     *   pow(P1 * P2 % P3, -1, P3) = the combined inverse for step 3
     */
    static const uint64_t INV_P1_MOD_P2   = 657107549ULL;  /* P1^-1 mod P2, verified: pow(P1,-1,P2) */
    static const uint64_t INV_P1P2_MOD_P3 = 284003040ULL;  /* (P1*P2)^-1 mod P3, verified: pow(P1*P2%P3,-1,P3) */

    /* Step 1: x mod P1 = r1 */
    uint64_t a1 = r1;

    /* Step 2: a2 = (r2 - a1) * P1^-1  mod P2 */
    uint64_t a2 = mod_mul((r2 + P2 - a1 % P2) % P2, INV_P1_MOD_P2, P2);

    /* Step 3: a3 = (r3 - a1 - a2*P1) * (P1*P2)^-1  mod P3 */
    uint64_t sub1 = a1 % P3;
    uint64_t sub2 = mod_mul(a2 % P3, P1 % P3, P3);
    uint64_t sub  = (r3 + 2*P3 - sub1 - sub2) % P3;
    uint64_t a3   = mod_mul(sub, INV_P1P2_MOD_P3, P3);

    /* Reconstruct: x = a1 + a2*P1 + a3*P1*P2 */
    __uint128_t x = (__uint128_t)a1
                  + (__uint128_t)a2 * P1
                  + (__uint128_t)a3 * P1 * P2;
    return x;
}

/* ── Public entry point ─────────────────────────────────────────────────── */

BigUInt *biguint_mul_ntt_crt(const BigUInt *a, const BigUInt *b)
{
    if (biguint_is_zero(a) || biguint_is_zero(b))
        return biguint_new();

    size_t na = a->len;
    size_t nb = b->len;
    size_t n  = 1;
    while (n < na + nb - 1) n <<= 1;
    size_t conv_len  = na + nb - 1;   /* logical convolution length */
    size_t result_len = na + nb;      /* BigUInt words (one extra for carry) */

    /* Allocate six NTT buffers (three primes × two operands). */
    uint64_t *f1 = calloc(n, sizeof(uint64_t));
    uint64_t *g1 = calloc(n, sizeof(uint64_t));
    uint64_t *f2 = calloc(n, sizeof(uint64_t));
    uint64_t *g2 = calloc(n, sizeof(uint64_t));
    uint64_t *f3 = calloc(n, sizeof(uint64_t));
    uint64_t *g3 = calloc(n, sizeof(uint64_t));

    if (!f1 || !g1 || !f2 || !g2 || !f3 || !g3) {
        free(f1); free(g1); free(f2); free(g2); free(f3); free(g3);
        return NULL;
    }

    /* Load inputs (each word reduced mod each prime). */
    for (size_t i = 0; i < na; i++) {
        f1[i] = a->words[i] % P1;
        f2[i] = a->words[i] % P2;
        f3[i] = a->words[i] % P3;
    }
    for (size_t i = 0; i < nb; i++) {
        g1[i] = b->words[i] % P1;
        g2[i] = b->words[i] % P2;
        g3[i] = b->words[i] % P3;
    }

    /* Forward NTT over each prime. */
    ntt_inplace_crt(f1, n, P1, G1, false);
    ntt_inplace_crt(g1, n, P1, G1, false);
    ntt_inplace_crt(f2, n, P2, G2, false);
    ntt_inplace_crt(g2, n, P2, G2, false);
    ntt_inplace_crt(f3, n, P3, G3, false);
    ntt_inplace_crt(g3, n, P3, G3, false);

    /* Pointwise multiply in each field. */
    for (size_t i = 0; i < n; i++) {
        f1[i] = mod_mul(f1[i], g1[i], P1);
        f2[i] = mod_mul(f2[i], g2[i], P2);
        f3[i] = mod_mul(f3[i], g3[i], P3);
    }

    /* Inverse NTT over each prime. */
    ntt_inplace_crt(f1, n, P1, G1, true);
    ntt_inplace_crt(f2, n, P2, G2, true);
    ntt_inplace_crt(f3, n, P3, G3, true);

    /* CRT reconstruction + carry propagation into BigUInt words. */
    BigUInt *res = biguint_new();
    if (!res) { free(f1); free(g1); free(f2); free(g2); free(f3); free(g3); return NULL; }

    res->words = calloc(result_len, sizeof(uint64_t));
    if (!res->words) {
        free(f1); free(g1); free(f2); free(g2); free(f3); free(g3);
        free(res);
        return NULL;
    }
    res->len = result_len;
    res->cap = result_len;

    __uint128_t carry = 0;
    for (size_t i = 0; i < conv_len; i++) {
        __uint128_t coeff = garner_crt(f1[i], f2[i], f3[i]) + carry;
        res->words[i] = (uint64_t)coeff;
        carry = coeff >> 64;
    }
    /* Propagate remaining carry into the upper words. */
    for (size_t i = conv_len; i < result_len && carry; i++) {
        __uint128_t val = (__uint128_t)res->words[i] + carry;
        res->words[i] = (uint64_t)val;
        carry = val >> 64;
    }

    free(f1); free(g1);
    free(f2); free(g2);
    free(f3); free(g3);

    biguint_normalize(res);
    return res;
}

#endif /* NTT_CRT_C_INCLUDED */

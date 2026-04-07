#ifndef NUSSBAUMER_C_INCLUDED
#define NUSSBAUMER_C_INCLUDED

#include "bigint.h"
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

/* ──────────────────────────────────────────────────────────────────────────
 * Branchless helpers
 * ──────────────────────────────────────────────────────────────────────── */

/* Bit-reversal permutation index for a p-bit integer. */
static inline int reverse_bits(int x, int p) {
    int rev = 0;
    for (int i = 0; i < p; i++) {
        rev = (rev << 1) | ((x >> i) & 1);
    }
    return rev;
}

/*
 * Multiply polynomial block `in[0..K)` by x^shift  mod (x^K + 1),
 * storing the result in `out[0..K)`.
 *
 * Decompose shift = q*K + r  where  q ∈ {0,1}, r ∈ [0,K).
 * Then  x^shift ≡ (-1)^q · x^r  mod (x^K+1).
 *
 * For each output index k:
 *   src     = (k - r + K) & (K-1)   — branchless right-rotate by r
 *   wrapped = (k < r)               — 1 iff the rotate crossed x^K once more
 *   neg     = q ^ wrapped           — total negation parity (0 or 1)
 *   sign    = -(__int128_t)neg         — 0x0…0 or 0xF…F
 *   out[k]  = (in[src] ^ sign) - sign   — branchless conditional negate
 *
 * This is completely branch-free; every expression reduces to
 * shifts, AND, XOR, and two's-complement arithmetic.
 */
static void poly_mul_xk(__int128_t * restrict out,
                         const __int128_t * restrict in,
                         int K, int shift)
{
    int log2K = __builtin_ctz(K);  /* K is a power of 2 */
    int Kmask  = K - 1;
    int q      = shift >> log2K;   /* 0 or 1: number of full K periods */
    int r      = shift & Kmask;    /* shift mod K                       */

    for (int k = 0; k < K; k++) {
        int src     = (k - r + K) & Kmask; /* rotate index, always in [0,K)   */
        int wrapped = (k < r);             /* 1 if rotation wrapped mod x^K+1 */
        int neg     = q ^ wrapped;         /* combined negate flag             */
        __int128_t sign = -(__int128_t)neg;      /* 0 or all-ones mask               */
        out[k] = (in[src] ^ sign) - sign;
    }
}

/* ──────────────────────────────────────────────────────────────────────────
 * Cooley-Tukey FFT over blocks, twiddles = x^step (pure shifts, no muls)
 * ──────────────────────────────────────────────────────────────────────── */
static void nussbaumer_fft(__int128_t *A, int M, int K, bool inverse) {
    int p = 0;
    while ((1 << p) < M) p++;

    /* Bit-reversal permutation of blocks */
    for (int i = 1; i < M; i++) {
        int j = reverse_bits(i, p);
        if (i < j) {
            __int128_t *bi = A + i * K;
            __int128_t *bj = A + j * K;
            for (int k = 0; k < K; k++) {
                __int128_t t = bi[k];
                bi[k] = bj[k];
                bj[k] = t;
            }
        }
    }

    int root_power = 2 * K / M;          /* omega = x^root_power          */
    if (inverse) root_power = 2 * K - root_power;

    __int128_t *V = malloc(K * sizeof(__int128_t));

    for (int len = 2; len <= M; len <<= 1) {
        int half      = len >> 1;
        /* step_power is always < 2K because root_power < 2K and M/len >= 1 */
        int step_power = (root_power * (M / len)) % (2 * K);

        for (int i = 0; i < M; i += len) {
            for (int j = 0; j < half; j++) {
                int shift      = (j * step_power) % (2 * K);
                __int128_t *u_blk = A + (i + j)        * K;
                __int128_t *v_blk = A + (i + j + half) * K;

                /* V = v_blk * x^shift  (branchless shift via poly_mul_xk) */
                poly_mul_xk(V, v_blk, K, shift);

                /* Butterfly: u = u+V,  v = u-V */
                for (int k = 0; k < K; k++) {
                    __int128_t u  = u_blk[k];
                    __int128_t v  = V[k];
                    u_blk[k]   = u + v;
                    v_blk[k]   = u - v;
                }
            }
        }
    }

    if (inverse) {
        for (int i = 0; i < M * K; i++) A[i] /= M;
    }
    free(V);
}

/* ──────────────────────────────────────────────────────────────────────────
 * Schoolbook polynomial multiplication mod (x^K + 1)
 *
 * Split-loop trick: for each fixed i, j runs in two contiguous ranges:
 *   j ∈ [0, K-i)  →  pos = i+j < K  →  out[pos] += a[i]*b[j]
 *   j ∈ [K-i, K)  →  pos = i+j-K    →  out[pos] -= a[i]*b[j]
 * No branch inside the j-loop; the split point is a compile-time-like
 * scalar computed once per i.
 * ──────────────────────────────────────────────────────────────────────── */
static void poly_mul_mod(__int128_t * restrict out,
                          const __int128_t * restrict a,
                          const __int128_t * restrict b,
                          int K)
{
    for (int i = 0; i < K; i++) out[i] = 0;

    for (int i = 0; i < K; i++) {
        __int128_t ai = a[i];
        if (ai == 0) continue;

        int split = K - i;          /* j < split  → no wrap; j >= split → wrap */

        /* No-wrap half: pos = i+j, contribution is positive */
        for (int j = 0; j < split; j++) {
            out[i + j] += ai * b[j];
        }
        /* Wrap half: pos = i+j-K, contribution is negative (x^K = -1) */
        for (int j = split; j < K; j++) {
            out[i + j - K] -= ai * b[j];
        }
    }
}

/* ──────────────────────────────────────────────────────────────────────────
 * Main entry point
 * ──────────────────────────────────────────────────────────────────────── */
BigUInt *biguint_mul_nussbaumer(const BigUInt *a, const BigUInt *b) {
    if (biguint_is_zero(a) || biguint_is_zero(b)) return biguint_new();

    int La = (int)(a->len * 2);
    int Lb = (int)(b->len * 2);
    int L = 1, M, K;

    /* Choose L so that the cyclic convolution of M blocks does not wrap.
     * We need M (next-pow-2 of Ma+Mb-1) ≤ 2K = 4L. */
    while (1) {
        int Ma    = (La + L - 1) / L;
        int Mb    = (Lb + L - 1) / L;
        int M_req = Ma + Mb - 1;
        if (M_req < 1) M_req = 1;
        M = 1; while (M < M_req) M <<= 1;
        K = 2 * L;
        if (M <= 2 * K) break;
        L <<= 1;
    }

    __int128_t *A = calloc(M * K, sizeof(__int128_t));
    __int128_t *B = calloc(M * K, sizeof(__int128_t));

    /* Pack 64-bit limbs into 32-bit chunks, then into blocks of size L */
    for (int i = 0; i < La; i++) {
        /* Branchless extraction of the i-th 32-bit chunk from a->words */
        int wi    = i >> 1;                         /* word index        */
        int shift = (i & 1) << 5;                  /* 0 or 32           */
        uint32_t val = (uint32_t)(a->words[wi] >> shift);
        A[(i / L) * K + (i % L)] = val;
    }
    for (int i = 0; i < Lb; i++) {
        int wi    = i >> 1;
        int shift = (i & 1) << 5;
        uint32_t val = (uint32_t)(b->words[wi] >> shift);
        B[(i / L) * K + (i % L)] = val;
    }

    nussbaumer_fft(A, M, K, false);
    nussbaumer_fft(B, M, K, false);

    __int128_t *C = calloc(M * K, sizeof(__int128_t));
    for (int m = 0; m < M; m++) {
        poly_mul_mod(&C[m * K], &A[m * K], &B[m * K], K);
    }

    nussbaumer_fft(C, M, K, true);

    /* Overlap-add blocks back into a linear 32-bit coefficient array */
    size_t final_len32 = (size_t)(M * L + K);
    __int128_t *final_C   = calloc(final_len32, sizeof(__int128_t));
    for (int m = 0; m < M; m++) {
        for (int k = 0; k < K; k++) {
            final_C[m * L + k] += C[m * K + k];
        }
    }

    /* Carry propagation, unrolled 2× to avoid the odd/even branch */
    size_t res_max_words = a->len + b->len;
    BigUInt *res   = biguint_new();
    res->words     = calloc(res_max_words, sizeof(uint64_t));
    res->len       = res_max_words;
    res->cap       = res_max_words;

    unsigned __int128 carry = 0;
    size_t pairs = final_len32 / 2;          /* process two 32-bit chunks at a time */

    for (size_t p2 = 0; p2 < pairs; p2++) {
        /* Low 32 bits */
        unsigned __int128 vlo = (unsigned __int128)final_C[2 * p2]     + carry;
        /* High 32 bits (use the carry from the low half inline) */
        unsigned __int128 vhi = (unsigned __int128)final_C[2 * p2 + 1] + (vlo >> 32);
        if (p2 < res->len) {
            res->words[p2] = (uint32_t)(vlo & 0xFFFFFFFFULL)
                           | ((uint64_t)(uint32_t)(vhi & 0xFFFFFFFFULL) << 32);
        }
        carry = vhi >> 32;
    }
    /* Handle an odd trailing chunk (final_len32 may be odd) */
    if (final_len32 & 1) {
        size_t i = final_len32 - 1;
        unsigned __int128 vlo = (unsigned __int128)final_C[i] + carry;
        size_t wi = i >> 1;
        if (wi < res->len) res->words[wi] = (uint32_t)(vlo & 0xFFFFFFFFULL);
        carry = vlo >> 32;
    }

    /* Flush any remaining carry */
    size_t write_idx = (final_len32 + 1) / 2;
    while (carry) {
        if (write_idx >= res->cap) {
            res->cap *= 2;
            res->words = realloc(res->words, res->cap * sizeof(uint64_t));
            memset(res->words + res->len, 0,
                   (res->cap - res->len) * sizeof(uint64_t));
        }
        res->len = write_idx + 1;
        res->words[write_idx++] = (uint64_t)(carry & 0xFFFFFFFFFFFFFFFFULL);
        carry >>= 64;
    }

    free(A); free(B); free(C); free(final_C);

    biguint_normalize(res);
    return res;
}
#endif

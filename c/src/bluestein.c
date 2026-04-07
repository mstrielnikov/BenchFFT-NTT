/**
 * bluestein.c — Bluestein / Chirp-Z Transform  for BigUInt multiplication
 *
 * The Chirp-Z Transform (CZT) computes a DFT of ARBITRARY length N using
 * an inner FFT of length M = next_pow2(2N−1).  This removes the standard
 * requirement that N must divide p−1 (NTT) or be a power of 2 (radix-2
 * FFT), and gives exact control over the evaluation length.
 *
 * Derivation — identity:  kn = (k² + n² − (k−n)²) / 2
 *
 *   X[k] = Σ x[n]·ω^{kn}          ω = exp(+2πi/N)
 *         = exp(+πi k²/N) · Σ [x[n]·exp(+πi n²/N)] · exp(−πi (k−n)²/N)
 *         = W_k  ·  (A ⊛ B)[k]
 *
 * where A[n] = x[n]·exp(+πi n²/N),  B[n] = exp(−πi n²/N)  (chirp),
 * and ⊛ is linear convolution computed via the inner power-of-2 FFT.
 *
 * Precision:
 *   Chunk size = 15 bits (max value 32767).
 *   Max convolution coefficient = N × 32767² ≤ 16384 × 1.07×10⁹ ≈ 1.75×10¹³
 *   IEEE-754 double mantissa holds exact integers up to 2^53 ≈ 9×10¹⁵.
 *   So 15-bit chunks are safe for N ≤ 16384 (≈ 245 760-bit inputs).
 *
 * Complexity vs. fft_split:
 *   For equal-size n-chunk inputs the inner FFT length is
 *     M_bluestein = next_pow2(4n−1) = 4n  (if n is a power of 2)
 *   vs. M_standard = 2n.  Bluestein therefore uses ~2× larger inner FFT.
 *   Its primary advantage is eliminating zero-padding when n is NOT a
 *   power of 2; for this benchmark (powers-of-2 word counts) performance
 *   is comparable to fft_split, demonstrating the algorithm.
 */

#ifndef BLUESTEIN_C_INCLUDED
#define BLUESTEIN_C_INCLUDED

#include "bigint.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ── Chunk size ─────────────────────────────────────────────────────────────
 * 15-bit chunks: each uint64_t word → 4 chunks (60 bits used, 4 wasted)   */
#define CZT_BITS  15
#define CZT_BASE  (1LL << CZT_BITS)   /* 32768 */

/* ═══════════════════════════════════════════════════════════════════════════
 * Internal iterative Cooley-Tukey FFT (split real/imaginary arrays)
 *
 * Convention (matching fft_split.c):
 *   inverse=false  →  forward  DFT: twiddle = exp(+2πi k/len)
 *   inverse=true   →  inverse IDFT: twiddle = exp(−2πi k/len), divide by n
 * ═══════════════════════════════════════════════════════════════════════════ */
static void czt_fft(double *re, double *im, size_t n, bool inverse)
{
    /* Bit-reversal permutation */
    for (size_t i = 1, j = 0; i < n; i++) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            double t;
            t = re[i]; re[i] = re[j]; re[j] = t;
            t = im[i]; im[i] = im[j]; im[j] = t;
        }
    }
    /* Butterfly stages */
    for (size_t len = 2; len <= n; len <<= 1) {
        double angle = (inverse ? -2.0 : +2.0) * M_PI / (double)len;
        double wl_re = cos(angle), wl_im = sin(angle);

        for (size_t i = 0; i < n; i += len) {
            double w_re = 1.0, w_im = 0.0;
            for (size_t k = 0; k < len / 2; k++) {
                size_t u = i + k, v = i + k + len / 2;
                double u_re = re[u], u_im = im[u];
                double vw_re = re[v] * w_re - im[v] * w_im;
                double vw_im = re[v] * w_im + im[v] * w_re;
                re[u] = u_re + vw_re;  im[u] = u_im + vw_im;
                re[v] = u_re - vw_re;  im[v] = u_im - vw_im;
                double nw_re = w_re * wl_re - w_im * wl_im;
                double nw_im = w_re * wl_im + w_im * wl_re;
                w_re = nw_re; w_im = nw_im;
            }
        }
    }
    if (inverse) {
        for (size_t i = 0; i < n; i++) { re[i] /= (double)n; im[i] /= (double)n; }
    }
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Bluestein / Chirp-Z DFT of arbitrary length N
 *
 *  x_re[0..N-1], x_im[0..N-1]  →  y_re[0..N-1], y_im[0..N-1]
 *
 * sign_chirp = +1 for forward, −1 for inverse.
 * Inverse additionally divides by N.
 * ═══════════════════════════════════════════════════════════════════════════ */
static void czt_dft(double * restrict y_re, double * restrict y_im,
                    const double * restrict x_re, const double * restrict x_im,
                    size_t N, double sign_chirp)
{
    /* Inner FFT length M ≥ 2N−1 (next power of 2) */
    size_t M = 1;
    while (M < 2 * N - 1) M <<= 1;

    double *ar = calloc(M, sizeof(double));  /* A_re */
    double *ai = calloc(M, sizeof(double));  /* A_im */
    double *br = calloc(M, sizeof(double));  /* B_re (chirp)  */
    double *bi = calloc(M, sizeof(double));  /* B_im          */

    /* Precompute chirp table and weight input:
     *   A[n] = x[n] · exp(+sign·πi·n²/N)
     *   B[n] = exp(−sign·πi·n²/N)          */
    for (size_t n = 0; n < N; n++) {
        double ang = sign_chirp * M_PI * (double)(n * n) / (double)N;
        double c = cos(ang), s = sin(ang);
        /* A[n] = x[n] · exp(sign·πi·n²/N) */
        ar[n] = x_re[n] * c - x_im[n] * s;
        ai[n] = x_re[n] * s + x_im[n] * c;
        /* B[n] = conj of A's chirp */
        br[n] =  c;
        bi[n] = -s;
    }
    /* B for "negative" indices (B is even): B[M−n] = B[n], n = 1..N−1 */
    for (size_t n = 1; n < N; n++) {
        br[M - n] = br[n];
        bi[M - n] = bi[n];
    }

    /* Inner FFT of A and B (forward: positive exponent) */
    czt_fft(ar, ai, M, false);
    czt_fft(br, bi, M, false);

    /* Pointwise multiply: C = A_fft · B_fft */
    for (size_t k = 0; k < M; k++) {
        double re = ar[k] * br[k] - ai[k] * bi[k];
        double im = ar[k] * bi[k] + ai[k] * br[k];
        ar[k] = re; ai[k] = im;
    }

    /* Inner IFFT to get convolution */
    czt_fft(ar, ai, M, true);

    /* Apply output chirp weight: y[k] = exp(sign·πi·k²/N) · conv[k] */
    for (size_t k = 0; k < N; k++) {
        double ang = sign_chirp * M_PI * (double)(k * k) / (double)N;
        double c = cos(ang), s = sin(ang);
        y_re[k] = ar[k] * c - ai[k] * s;
        y_im[k] = ar[k] * s + ai[k] * c;
    }

    /* Inverse: divide by N */
    if (sign_chirp < 0) {
        for (size_t k = 0; k < N; k++) {
            y_re[k] /= (double)N;
            y_im[k] /= (double)N;
        }
    }

    free(ar); free(ai); free(br); free(bi);
}

/* ═══════════════════════════════════════════════════════════════════════════
 * BigUInt multiplication via Bluestein / Chirp-Z
 *
 * Steps:
 *   1. Split each uint64_t word into ⌊64/CZT_BITS⌋ = 4 chunks (15-bit each).
 *   2. Zero-pad both chunk arrays to N = na + nb (chosen as the exact DFT
 *      length for linear convolution — no extra power-of-2 rounding).
 *   3. Forward CZT length-N DFT of both.
 *   4. Pointwise complex multiply.
 *   5. Inverse CZT length-N DFT.
 *   6. Round real parts to int64_t; carry-propagate in base CZT_BASE.
 * ═══════════════════════════════════════════════════════════════════════════ */
BigUInt *biguint_mul_bluestein(const BigUInt *a, const BigUInt *b)
{
    if (biguint_is_zero(a) || biguint_is_zero(b)) return biguint_new();

    /* --- 1. Split into 15-bit chunks --------------------------------------- */
    /* Number of chunks = ceil(bits / CZT_BITS) = ceil(len*64 / CZT_BITS) */
    size_t na = (a->len * 64 + CZT_BITS - 1) / CZT_BITS;
    size_t nb = (b->len * 64 + CZT_BITS - 1) / CZT_BITS;
    size_t N  = na + nb;   /* exact DFT length for linear convolution */

    double *xre = calloc(N, sizeof(double));
    double *xim = calloc(N, sizeof(double));
    double *yre = calloc(N, sizeof(double));
    double *yim = calloc(N, sizeof(double));
    double *Are = calloc(N, sizeof(double));
    double *Aim = calloc(N, sizeof(double));
    double *Bre = calloc(N, sizeof(double));
    double *Bim = calloc(N, sizeof(double));
    double *Cre = calloc(N, sizeof(double));
    double *Cim = calloc(N, sizeof(double));

    /* Extract chunks for A */
    for (size_t k = 0; k < na; k++) {
        size_t bit_lo = k * CZT_BITS;
        size_t wi = bit_lo / 64;
        size_t bi_off = bit_lo % 64;
        uint64_t val = (wi < a->len) ? (a->words[wi] >> bi_off) : 0;
        if (bi_off + CZT_BITS > 64 && wi + 1 < a->len)
            val |= a->words[wi + 1] << (64 - bi_off);
        xre[k] = (double)(val & (CZT_BASE - 1));
    }
    /* Extract chunks for B */
    for (size_t k = 0; k < nb; k++) {
        size_t bit_lo = k * CZT_BITS;
        size_t wi = bit_lo / 64;
        size_t bi_off = bit_lo % 64;
        uint64_t val = (wi < b->len) ? (b->words[wi] >> bi_off) : 0;
        if (bi_off + CZT_BITS > 64 && wi + 1 < b->len)
            val |= b->words[wi + 1] << (64 - bi_off);
        yre[k] = (double)(val & (CZT_BASE - 1));
    }

    /* --- 3. Forward CZT of both ------------------------------------------- */
    czt_dft(Are, Aim, xre, xim, N, +1.0);   /* forward: sign = +1 */
    czt_dft(Bre, Bim, yre, yim, N, +1.0);

    /* --- 4. Pointwise multiply -------------------------------------------- */
    for (size_t k = 0; k < N; k++) {
        Cre[k] = Are[k] * Bre[k] - Aim[k] * Bim[k];
        Cim[k] = Are[k] * Bim[k] + Aim[k] * Bre[k];
    }

    /* --- 5. Inverse CZT --------------------------------------------------- */
    czt_dft(xre, xim, Cre, Cim, N, -1.0);   /* inverse: sign = -1 */

    /* --- 6. Round to int64 and carry-propagate in base CZT_BASE ------------ */
    /* Max output length: ceil((na + nb) * CZT_BITS / 64) words */
    size_t rwords = (N * CZT_BITS + 63) / 64 + 2;
    uint64_t *r = calloc(rwords, sizeof(uint64_t));

    /* Accumulate into a big integer.
     * Each coefficient xre[k] ≈ C_k (an integer in base-CZT_BASE).
     * We pack CZT_BITS bits per coefficient into the output word array.
     * Use a 64-bit accumulator with carry to avoid double-precision subtleties. */
    uint64_t carry = 0;
    for (size_t k = 0; k < N; k++) {
        int64_t coeff = (int64_t)round(xre[k]) + (int64_t)carry;
        carry = 0;
        if (coeff < 0) {
            /* Should not happen for valid (non-negative) BigUInt products */
            coeff = 0;
        }
        /* Extract CZT_BITS of this coefficient into the bit stream */
        carry     = (uint64_t)coeff >> CZT_BITS;
        uint64_t digit = (uint64_t)coeff & (CZT_BASE - 1);

        /* Place digit into the output bit array at bit position k*CZT_BITS */
        size_t bit_lo = k * CZT_BITS;
        size_t wi     = bit_lo / 64;
        size_t bi_off = bit_lo % 64;
        if (wi < rwords)         r[wi]     |= digit << bi_off;
        if (bi_off + CZT_BITS > 64 && wi + 1 < rwords)
            r[wi + 1] |= digit >> (64 - bi_off);
    }
    /* Flush remaining carry */
    {
        size_t bit_lo = N * CZT_BITS;
        size_t wi     = bit_lo / 64;
        size_t bi_off = bit_lo % 64;
        if (carry && wi < rwords)       r[wi]     |= carry << bi_off;
        if (carry && bi_off + 64 > 64 && wi + 1 < rwords)
            r[wi + 1] |= carry >> (64 - bi_off);
    }

    free(xre); free(xim); free(yre); free(yim);
    free(Are); free(Aim); free(Bre); free(Bim);
    free(Cre); free(Cim);

    BigUInt *res = biguint_new();
    res->words = r;
    res->len   = rwords;
    res->cap   = rwords;
    biguint_normalize(res);
    return res;
}

#endif /* BLUESTEIN_C_INCLUDED */

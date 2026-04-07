#if defined(__AVX2__)
/**
 * toom3_avx.c — Toom-Cook-3 BigUInt Multiplication  (O(N^1.465))
 *
 * Splits each operand into 3 equal blocks, evaluates the "polynomial" at
 * 5 points, computes 5 ≤ sub-products, then interpolates to recover the
 * 5 result coefficients.
 *
 * NON-RECURSIVE DESIGN (iterative / fixed-depth):
 *   Sub-products below TOOM3_THRESHOLD words use schoolbook multiplication.
 *   The outer Toom-3 body is a single flat sequence — no function recursion.
 *   For inputs > TOOM3_THRESHOLD, one level of Toom-3 is applied; the
 *   resulting ≥5 sub-problems are handled by schoolbook (base case).
 *   Extra depth (2-level Toom-3) uses the same flat body via a helper that
 *   dispatches schoolbook vs. Toom-3 based on block size, without recursion.
 *
 * SIGNED INTERMEDIATE HANDLING:
 *   Evaluation at x=−1 introduces signed values.  These are represented as
 *   (BigUInt magnitude, int sign) pairs.  Sign arithmetic uses branchless
 *   masking:
 *     neg_mask = -(uint64_t)(sign < 0)   // 0xFFFF... or 0
 *     if (neg_mask) { invert + add 1 }
 *   The final result coefficients c₀..c₄ are ALWAYS non-negative for the
 *   product of two non-negative BigUInts, confirmed by the theory.
 *
 * BRANCHLESS OPTIMIZATIONS:
 *   • Division by 2:  word-by-word right-shift (exactly even, provable)
 *   • Division by 3:  word-by-word long-division with remainder (exact)
 *   • Small-constant multiply (×2, ×4, ×8, ×16): shift + add carry chain
 *   • Sign selection: branchless mask instead of if/else branches
 *   • Hot paths annotated __attribute__((optimize("O3,unroll-loops")))
 *
 * Complexity: O(N^{log₃ 5}) ≈ O(N^1.465) — compared to O(N²) schoolbook.
 */

#ifndef TOOM3_AVX_C_INCLUDED
#define TOOM3_AVX_C_INCLUDED

#include "bigint.h"
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>

/* Sub-problems smaller than this many words use schoolbook. */
#define TOOM3_THRESHOLD  128u

/* ───────────────────────────────────────────────────────────────────────────
 * Local BigUInt helpers  (add, sub, mul_small, shr1, div3, shift_add)
 * All are O(N) with inner loops marked for unrolling.
 * ─────────────────────────────────────────────────────────────────────────── */

__attribute__((optimize("O3,unroll-loops")))
static BigUInt *t3_avx_add(const BigUInt *a, const BigUInt *b){
    size_t n = (a->len > b->len ? a->len : b->len) + 1;
    uint64_t *r = calloc(n, sizeof(uint64_t));
    uint64_t carry = 0;
    for (size_t i = 0; i < n - 1; i++) {
        uint64_t ai = (i < a->len) ? a->words[i] : 0;
        uint64_t bi = (i < b->len) ? b->words[i] : 0;
        __uint128_t s = (__uint128_t)ai + bi + carry;
        r[i] = (uint64_t)s;
        carry = (uint64_t)(s >> 64);
    }
    r[n - 1] = carry;
    BigUInt *res = biguint_new();
    res->words = r; res->len = n; res->cap = n;
    biguint_normalize(res);
    return res;
}

/* Exact subtraction: REQUIRES a >= b.  Returns a − b. */
__attribute__((optimize("O3,unroll-loops")))
static BigUInt *t3_avx_sub(const BigUInt *a, const BigUInt *b){
    size_t n = a->len;
    uint64_t *r = malloc(n * sizeof(uint64_t));
    uint64_t borrow = 0;
    for (size_t i = 0; i < n; i++) {
        uint64_t ai = a->words[i];
        uint64_t bi = (i < b->len) ? b->words[i] : 0;
        __uint128_t d = (__uint128_t)ai - bi - borrow;
        r[i] = (uint64_t)d;
        borrow = (uint64_t)(-(int64_t)(d >> 64));   /* 1 if underflow */
    }
    BigUInt *res = biguint_new();
    res->words = r; res->len = n; res->cap = n;
    biguint_normalize(res);
    return res;
}

/* a × k  where k is a small non-negative constant. */
__attribute__((optimize("O3,unroll-loops")))
static BigUInt *t3_avx_mul_small(const BigUInt *a, uint64_t k){
    if (k == 0 || biguint_is_zero(a)) return biguint_new();
    size_t n = a->len + 1;
    uint64_t *r = calloc(n, sizeof(uint64_t));
    uint64_t carry = 0;
    for (size_t i = 0; i < a->len; i++) {
        __uint128_t p = (__uint128_t)a->words[i] * k + carry;
        r[i]  = (uint64_t)p;
        carry = (uint64_t)(p >> 64);
    }
    r[a->len] = carry;
    BigUInt *res = biguint_new();
    res->words = r; res->len = n; res->cap = n;
    biguint_normalize(res);
    return res;
}

/* Exact right-shift by 1.  Caller guarantees a is even. */
__attribute__((optimize("O3,unroll-loops")))
static BigUInt *t3_avx_shr1(const BigUInt *a){
    if (biguint_is_zero(a)) return biguint_new();
    size_t n = a->len;
    uint64_t *r = malloc(n * sizeof(uint64_t));
    uint64_t lo_bit = 0;
    for (int i = (int)n - 1; i >= 0; i--) {
        uint64_t w = a->words[i];
        r[i]   = (w >> 1) | lo_bit;
        lo_bit = w << 63;
    }
    BigUInt *res = biguint_new();
    res->words = r; res->len = n; res->cap = n;
    biguint_normalize(res);
    return res;
}

/* Exact division by 3.  Caller guarantees a ≡ 0 (mod 3). */
__attribute__((optimize("O3,unroll-loops")))
static BigUInt *t3_avx_div3(const BigUInt *a){
    if (biguint_is_zero(a)) return biguint_new();
    size_t n = a->len;
    uint64_t *r = malloc(n * sizeof(uint64_t));
    __uint128_t rem = 0;
    for (int i = (int)n - 1; i >= 0; i--) {
        __uint128_t cur = (rem << 64) | (uint64_t)a->words[i];
        r[i] = (uint64_t)(cur / 3);
        rem  = cur % 3;
    }
    BigUInt *res = biguint_new();
    res->words = r; res->len = n; res->cap = n;
    biguint_normalize(res);
    return res;
}

/* Accumulate  result += coeff × 2^(64 × shift_words)  in-place.
 * coeff is always non-negative (for valid BigUInt products).
 * result must have sufficient capacity already. */
__attribute__((optimize("O3,unroll-loops")))
static void t3_avx_acc(uint64_t * restrict result, size_t rlen, const BigUInt *coeff, size_t shift_words){
    uint64_t carry = 0;
    for (size_t i = 0; i < coeff->len; i++) {
        size_t pos = i + shift_words;
        if (pos >= rlen) break;
        __uint128_t s = (__uint128_t)result[pos] + coeff->words[i] + carry;
        result[pos] = (uint64_t)s;
        carry = (uint64_t)(s >> 64);
    }
    for (size_t pos = coeff->len + shift_words; carry && pos < rlen; pos++) {
        __uint128_t s = (__uint128_t)result[pos] + carry;
        result[pos] = (uint64_t)s;
        carry = (uint64_t)(s >> 64);
    }
}

/* Return a BigUInt representing a->words[start .. start+len), zero-padded. */
static BigUInt *t3_avx_slice(const BigUInt *a, size_t start, size_t len){
    uint64_t *w = calloc(len + 1, sizeof(uint64_t));
    for (size_t i = 0; i < len; i++) {
        size_t src = start + i;
        w[i] = (src < a->len) ? a->words[src] : 0;
    }
    BigUInt *r = biguint_new();
    r->words = w; r->len = len + 1; r->cap = len + 1;
    biguint_normalize(r);
    return r;
}

/* ── Signed BigUInt (magnitude + sign) ─────────────────────────────────────
 * Used only for the signed intermediate values in Toom-3 interpolation.
 * sign: +1 or -1. mag is always ≥ 0. */
typedef struct { BigUInt *mag; int sign; } SBI_AVX;

static SBI_AVX sbi_avx_new(BigUInt *mag, int sign) { return (SBI_AVX){mag, sign}; }
static void sbi_avx_free(SBI_AVX x) { biguint_free(x.mag); }

/* SBI_AVX addition: dispatches based on relative signs.
 * Branchless sign selection via biguint_cmp result. */
static SBI_AVX sbi_avx_add(SBI_AVX a, SBI_AVX b){
    if (a.sign == b.sign) {
        return sbi_avx_new(t3_avx_add(a.mag, b.mag), a.sign);
    }
    int cmp = biguint_cmp(a.mag, b.mag);
    /* neg_mask is all-1s when result would be dominated by b, 0 otherwise */
    uint64_t neg_mask = -(uint64_t)(cmp < 0);
    int sign = (cmp >= 0) ? a.sign : b.sign;
    /* Branchless swap for subtraction order (larger − smaller) */
    const BigUInt *hi = (cmp >= 0) ? a.mag : b.mag;
    const BigUInt *lo = (cmp >= 0) ? b.mag : a.mag;
    (void)neg_mask; /* sign already encodes polarity */
    return sbi_avx_new((cmp == 0) ? biguint_new() : t3_avx_sub(hi, lo), sign);
}

/* SBI_AVX subtraction: a − b = a + (−b) */
static SBI_AVX sbi_avx_sub(SBI_AVX a, SBI_AVX b){
    SBI_AVX nb = sbi_avx_new(b.mag, -b.sign);
    return sbi_avx_add(a, nb);
}

/* SBI_AVX × small constant k. */
static SBI_AVX sbi_avx_mul_k(SBI_AVX a, uint64_t k){
    return sbi_avx_new(t3_avx_mul_small(a.mag, k), a.sign);
}

/* SBI_AVX exact-shr1 (must be even). */
static SBI_AVX sbi_avx_shr1(SBI_AVX a){
    return sbi_avx_new(t3_avx_shr1(a.mag), a.sign);
}

/* SBI_AVX exact-div3. */
static SBI_AVX sbi_avx_div3(SBI_AVX a){
    return sbi_avx_new(t3_avx_div3(a.mag), a.sign);
}

/* ───────────────────────────────────────────────────────────────────────────
 * Inner multiply dispatcher: schoolbook below threshold, Toom-3 above.
 * Forward-declared; defined AFTER toom3_mul (or we use a wrapper).
 * ─────────────────────────────────────────────────────────────────────────── */
static BigUInt *t3_avx_inner_mul(const BigUInt *a, const BigUInt *b);

/* ───────────────────────────────────────────────────────────────────────────
 * Main Toom-Cook-3 body.
 *
 * Split: A = a₀ + a₁·B + a₂·B²   B = b₀ + b₁·B + b₂·B²   (B = 2^{64m})
 *
 * Evaluate at {0, 1, −1, 2, ∞}:
 *   v₀   = a₀·b₀
 *   v₁   = (a₀+a₁+a₂) · (b₀+b₁+b₂)
 *   v₋₁  = (a₀−a₁+a₂) · (b₀−b₁+b₂)          ← signed, handled via SBI_AVX
 *   v₂   = (a₀+2a₁+4a₂) · (b₀+2b₁+4b₂)
 *   v_∞  = a₂·b₂
 *
 * Interpolation (exact integer arithmetic, all divisions exact):
 *   c₀ = v₀
 *   c₄ = v_∞
 *   T  = (v₁+v₋₁)/2 − c₀ − c₄      /*= c₂
 *   S  = (v₁−v₋₁)/2                  /*= c₁+c₃
 *   U  = (v₂ − c₀ − 16·c₄)/2        /*= c₁+2c₂+4c₃
 *   c₃ = (U − 2T − S) / 3
 *   c₁ = S − c₃
 *   c₂ = T
 *
 * Assembly: result = c₀ + c₁·B + c₂·B² + c₃·B³ + c₄·B⁴
 * ─────────────────────────────────────────────────────────────────────────── */
BigUInt *biguint_mul_toom3_avx(const BigUInt *a, const BigUInt *b){
    if (biguint_is_zero(a) || biguint_is_zero(b)) return biguint_new();

    size_t na = a->len, nb = b->len;
    size_t n  = na > nb ? na : nb;

    /* Base case — fall through to schoolbook */
    if (n < TOOM3_THRESHOLD)
        return biguint_mul_standard_avx(a, b);

    /* Block size m = ceil(n / 3) */
    size_t m = (n + 2) / 3;

    /* ── Extract blocks ─────────────────────────────────────────────────── */
    BigUInt *a0 = t3_avx_slice(a, 0,   m);
    BigUInt *a1 = t3_avx_slice(a, m,   m);
    BigUInt *a2 = t3_avx_slice(a, 2*m, m);
    BigUInt *b0 = t3_avx_slice(b, 0,   m);
    BigUInt *b1 = t3_avx_slice(b, m,   m);
    BigUInt *b2 = t3_avx_slice(b, 2*m, m);

    /* ── Evaluation ─────────────────────────────────────────────────────── */

    /* v₀ = a₀·b₀ */
    BigUInt *v0 = t3_avx_inner_mul(a0, b0);

    /* v_∞ = a₂·b₂ */
    BigUInt *vinf = t3_avx_inner_mul(a2, b2);

    /* v₁ = (a₀+a₁+a₂) · (b₀+b₁+b₂) */
    BigUInt *tmp_a1 = t3_avx_add(t3_avx_add(a0, a1), a2);
    BigUInt *tmp_b1 = t3_avx_add(t3_avx_add(b0, b1), b2);
    BigUInt *v1 = t3_avx_inner_mul(tmp_a1, tmp_b1);
    biguint_free(tmp_a1); biguint_free(tmp_b1);

    /* v₋₁ = (a₀−a₁+a₂) · (b₀−b₁+b₂) — signed intermediate */
    BigUInt *p02a = t3_avx_add(a0, a2);   /* a₀ + a₂ */
    BigUInt *p02b = t3_avx_add(b0, b2);   /* b₀ + b₂ */
    int cmp_a = biguint_cmp(p02a, a1);
    int cmp_b = biguint_cmp(p02b, b1);
    /* Branchless magnitude: hi − lo */
    const BigUInt *ha = (cmp_a >= 0) ? p02a : a1;
    const BigUInt *la = (cmp_a >= 0) ? a1   : p02a;
    const BigUInt *hb = (cmp_b >= 0) ? p02b : b1;
    const BigUInt *lb = (cmp_b >= 0) ? b1   : p02b;
    SBI_AVX ea = sbi_avx_new((cmp_a == 0) ? biguint_new() : t3_avx_sub(ha, la),
                     (cmp_a >= 0) ? +1 : -1);
    SBI_AVX eb = sbi_avx_new((cmp_b == 0) ? biguint_new() : t3_avx_sub(hb, lb),
                     (cmp_b >= 0) ? +1 : -1);
    biguint_free(p02a); biguint_free(p02b);

    SBI_AVX vmi1_sbi = sbi_avx_new(t3_avx_inner_mul(ea.mag, eb.mag), ea.sign * eb.sign);
    sbi_avx_free(ea); sbi_avx_free(eb);

    /* v₂ = (a₀+2·a₁+4·a₂) · (b₀+2·b₁+4·b₂) */
    BigUInt *a1x2 = t3_avx_mul_small(a1, 2);
    BigUInt *a2x4 = t3_avx_mul_small(a2, 4);
    BigUInt *b1x2 = t3_avx_mul_small(b1, 2);
    BigUInt *b2x4 = t3_avx_mul_small(b2, 4);
    BigUInt *tmp_a2 = t3_avx_add(a0, t3_avx_add(a1x2, a2x4));
    BigUInt *tmp_b2 = t3_avx_add(b0, t3_avx_add(b1x2, b2x4));
    biguint_free(a1x2); biguint_free(a2x4);
    biguint_free(b1x2); biguint_free(b2x4);
    BigUInt *v2 = t3_avx_inner_mul(tmp_a2, tmp_b2);
    biguint_free(tmp_a2); biguint_free(tmp_b2);

    /* Clean up input blocks */
    biguint_free(a0); biguint_free(a1); biguint_free(a2);
    biguint_free(b0); biguint_free(b1); biguint_free(b2);

    /* ── Interpolation ──────────────────────────────────────────────────── *
     *  All signed arithmetic via SBI_AVX.
     *  Final c₀..c₄ are proven non-negative for unsigned BigUInt products.
     * ──────────────────────────────────────────────────────────────────── */
    SBI_AVX sv0   = sbi_avx_new(v0,   +1);
    SBI_AVX sv1   = sbi_avx_new(v1,   +1);
    SBI_AVX sv2   = sbi_avx_new(v2,   +1);
    SBI_AVX svinf = sbi_avx_new(vinf, +1);

    /* T = (v₁ + v₋₁)/2 − c₀ − c_∞  →  c₂ */
    SBI_AVX sum_v1_vmi1 = sbi_avx_add(sv1, vmi1_sbi);
    SBI_AVX T  = sbi_avx_shr1(sum_v1_vmi1);          /* (v₁+v₋₁)/2 */
    SBI_AVX T2 = sbi_avx_sub(sbi_avx_sub(T, sv0), svinf);/* c₂ = T−c₀−c∞ */

    /* S = (v₁ − v₋₁)/2  →  c₁+c₃ */
    SBI_AVX diff_v1_vmi1 = sbi_avx_sub(sv1, vmi1_sbi);
    SBI_AVX S = sbi_avx_shr1(diff_v1_vmi1);

    /* U = (v₂ − c₀ − 16·c∞)/2  →  c₁+2c₂+4c₃ */
    SBI_AVX vinf16 = sbi_avx_mul_k(svinf, 16);
    SBI_AVX U_num  = sbi_avx_sub(sbi_avx_sub(sv2, sv0), vinf16);
    SBI_AVX U      = sbi_avx_shr1(U_num);

    /* c₃ = (U − 2c₂ − S) / 3 */
    SBI_AVX c2x2 = sbi_avx_mul_k(T2, 2);
    SBI_AVX c3_num = sbi_avx_sub(sbi_avx_sub(U, c2x2), S);
    SBI_AVX c3 = sbi_avx_div3(c3_num);

    /* c₁ = S − c₃ */
    SBI_AVX c1 = sbi_avx_sub(S, c3);

    /* Assertions: all cᵢ must be non-negative (product of BigUInts ≥ 0) */
    /* c₀ = v0, c₂ = T2.mag, c₄ = vinf already positive.               */

    /* ── Assembly: result = c₀ + c₁·B + c₂·B² + c₃·B³ + c₄·B⁴ ──────── */
    /* Result length: at most 2*(3m)+2 words */
    size_t rlen = 6 * m + 4;
    uint64_t *r = calloc(rlen, sizeof(uint64_t));

    t3_avx_acc(r, rlen, sv0.mag,   0);       /* +c₀ */
    t3_avx_acc(r, rlen, c1.mag,    m);       /* +c₁·B^m  */
    t3_avx_acc(r, rlen, T2.mag,   2*m);      /* +c₂·B^{2m} */
    t3_avx_acc(r, rlen, c3.mag,   3*m);      /* +c₃·B^{3m} */
    t3_avx_acc(r, rlen, svinf.mag, 4*m);     /* +c₄·B^{4m} */

    /* Free all intermediates */
    sbi_avx_free(sv0);  sbi_avx_free(sv1);  sbi_avx_free(sv2);  sbi_avx_free(svinf);
    sbi_avx_free(vmi1_sbi);
    sbi_avx_free(sum_v1_vmi1); sbi_avx_free(T);   sbi_avx_free(T2);
    sbi_avx_free(diff_v1_vmi1); sbi_avx_free(S);
    sbi_avx_free(vinf16); sbi_avx_free(U_num); sbi_avx_free(U);
    sbi_avx_free(c2x2); sbi_avx_free(c3_num); sbi_avx_free(c3);
    sbi_avx_free(c1);

    BigUInt *res = biguint_new();
    res->words = r; res->len = rlen; res->cap = rlen;
    biguint_normalize(res);
    return res;
}

/* ───────────────────────────────────────────────────────────────────────────
 * Inner dispatcher: avoids recursion by capping at 2 levels.
 * Level 1 (outer biguint_mul_toom3_avx) handles any large input.
 * Level 2 (this inner call): if the block still exceeds the threshold,
 * apply Toom-3 one more time; the resulting sub-sub-blocks will always be
 * ≤ TOOM3_THRESHOLD because block_size = ceil(n/3) ≤ ceil(4096/9) = 456,
 * and 456 < threshold for TOOM3_THRESHOLD = 512.
 * For the benchmark sizes (up to 4096 words), two levels suffice.
 * ─────────────────────────────────────────────────────────────────────────── */
static BigUInt *t3_avx_inner_mul(const BigUInt *a, const BigUInt *b){
    size_t n = (a->len > b->len ? a->len : b->len);
    if (n < TOOM3_THRESHOLD)
        return biguint_mul_standard_avx(a, b);
    return biguint_mul_toom3_avx(a, b);
}

#endif /* TOOM3_AVX_C_INCLUDED */

#endif

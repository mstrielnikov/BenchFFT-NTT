/*
 * m61_math.h — Shared F_{M61} arithmetic primitives.
 *
 * Include-guarded so it is safe to pull in from multiple source files
 * in the same C translation unit (e.g. fft_mersenne.c + ntt_mersenne.c
 * both textually included into a single test binary).
 */
#pragma once

#include <stdint.h>

#define M61_MOD (2305843009213693951ULL)   /* 2^61 - 1 */
#define M61_P   (61)

static inline uint64_t m61_reduce(uint64_t x) {
    uint64_t t = (x >> M61_P) + (x & M61_MOD);
    t += (t >> M61_P);
    return t & M61_MOD;
}

static inline uint64_t m61_add(uint64_t a, uint64_t b) {
    uint64_t x = a + b;
    x += (x >> M61_P);
    return x & M61_MOD;
}

static inline uint64_t m61_sub(uint64_t a, uint64_t b) {
    uint64_t x = a + M61_MOD - b;
    x += (x >> M61_P);
    return x & M61_MOD;
}

static inline uint64_t m61_mul(uint64_t a, uint64_t b) {
    __uint128_t p  = (__uint128_t)a * b;
    uint64_t    hi = (uint64_t)(p >> M61_P);
    uint64_t    lo = (uint64_t)(p & M61_MOD);
    uint64_t    t  = hi + lo;
    t += (t >> M61_P);
    return t & M61_MOD;
}

static inline uint64_t m61_pow(uint64_t a, uint64_t e) {
    uint64_t res = 1;
    a &= M61_MOD;
    while (e) {
        if (e & 1) res = m61_mul(res, a);
        a = m61_mul(a, a);
        e >>= 1;
    }
    return res;
}

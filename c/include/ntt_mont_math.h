/*
 * ntt_mont_math.h — Shared generic NTT modular arithmetic.
 * Used by ntt_mont.c (scalar) and ntt_mont_avx.c (AVX2).
 */
#pragma once

#include <stdint.h>

static inline uint64_t mod_pow(uint64_t a, uint64_t e, uint64_t mod) {
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

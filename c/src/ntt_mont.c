#pragma once
#include <bigint.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <ntt_mont_math.h>


static void ntt_mont_inplace(uint64_t *a, size_t n, uint64_t mod, uint64_t root, bool inverse) {
    for (size_t i = 1, j = 0; i < n; i++) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            uint64_t t = a[i]; a[i] = a[j]; a[j] = t;
        }
    }
    
    uint64_t *roots = malloc(n * sizeof(uint64_t));
    roots[0] = 1;
    for (size_t len = 2; len <= n; len <<= 1) {
        uint64_t wlen = mod_pow(root, (mod - 1) / len, mod);
        if (inverse) wlen = mod_pow(wlen, mod - 2, mod);
        for (size_t i = len/2; i < len; i++) {
            roots[i] = mod_mul(roots[i - len/2], wlen, mod);
        }
    }
    
    size_t len = 2;
    while (len <= n) {
        for (size_t i = 0; i < n; i += len) {
            for (size_t j = 0; j < len / 2; j++) {
                uint64_t u = a[i + j];
                uint64_t v = mod_mul(a[i + j + len / 2], roots[j], mod);
                a[i + j] = u + v;
                if (a[i + j] >= mod) a[i + j] -= mod;
                a[i + j + len / 2] = u + mod - v;
                if (a[i + j + len / 2] >= mod) a[i + j + len / 2] -= mod;
            }
        }
        len <<= 1;
    }
    
    if (inverse) {
        uint64_t inv_n = mod_pow(n, mod - 2, mod);
        for (size_t i = 0; i < n; i++) {
            a[i] = mod_mul(a[i], inv_n, mod);
        }
    }
    
    free(roots);
}

BigUInt *biguint_mul_ntt_mont(const BigUInt *a, const BigUInt *b) {
    if (biguint_is_zero(a) || biguint_is_zero(b)) {
        return biguint_new();
    }
    
    const uint64_t MOD = 998244353;
    const uint64_t ROOT = 3;
    
    size_t n = 1;
    while (n < a->len + b->len - 1) n <<= 1;
    size_t result_len = a->len + b->len - 1;
    
    // Pre-allocate all vectors upfront (like Rust)
    uint64_t *fa = calloc(n, sizeof(uint64_t));
    uint64_t *fb = calloc(n, sizeof(uint64_t));
    
    if (!fa || !fb) {
        free(fa); free(fb);
        return NULL;
    }
    
    // Copy input words (modulo MOD)
    for (size_t i = 0; i < a->len; i++) fa[i] = a->words[i] % MOD;
    for (size_t i = 0; i < b->len; i++) fb[i] = b->words[i] % MOD;
    
    // NTT forward
    ntt_mont_inplace(fa, n, MOD, ROOT, false);
    ntt_mont_inplace(fb, n, MOD, ROOT, false);
    
    // Pointwise multiply
    for (size_t i = 0; i < n; i++) {
        fa[i] = mod_mul(fa[i], fb[i], MOD);
    }
    
    // Inverse NTT
    ntt_mont_inplace(fa, n, MOD, ROOT, true);
    
    // Pre-allocate result with exact size (like Rust Vec::with_capacity)
    BigUInt *res = biguint_new();
    if (!res) {
        free(fa); free(fb);
        return NULL;
    }
    
    // Pre-allocate result words array
    res->words = calloc(result_len, sizeof(uint64_t));
    if (!res->words) {
        free(fa); free(fb);
        free(res);
        return NULL;
    }
    
    // Copy result (only up to result_len)
    for (size_t i = 0; i < result_len; i++) {
        res->words[i] = fa[i];
    }
    res->len = result_len;
    res->cap = result_len;
    
    free(fa); free(fb);
    biguint_normalize(res);
    
    return res;
}

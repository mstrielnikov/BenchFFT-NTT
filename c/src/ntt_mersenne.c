#include <bigint.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#define M61_MOD (2305843009213693951ULL)
#define M61_P (61)

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

static uint64_t m61_mod_pow(uint64_t a, uint64_t e) {
    uint64_t res = 1;
    while (e) {
        if (e & 1) res = m61_reduce(res * a);
        a = m61_reduce(a * a);
        e >>= 1;
    }
    return res;
}

static inline uint64_t m61_mul(uint64_t a, uint64_t b) {
    return m61_reduce(a * b);
}

static void ntt_inplace_m61(uint64_t *a, size_t n, uint64_t root, bool inverse) {
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
        uint64_t wlen = m61_mod_pow(root, (M61_MOD - 1) / len);
        if (inverse) wlen = m61_mod_pow(wlen, M61_MOD - 2);
        for (size_t i = len/2; i < len; i++) {
            roots[i] = m61_mul(roots[i - len/2], wlen);
        }
    }
    
    size_t len = 2;
    while (len <= n) {
        for (size_t i = 0; i < n; i += len) {
            for (size_t j = 0; j < len / 2; j++) {
                uint64_t u = a[i + j];
                uint64_t v = m61_mul(a[i + j + len / 2], roots[j]);
                a[i + j] = m61_add(u, v);
                a[i + j + len / 2] = m61_sub(u, v);
            }
        }
        len <<= 1;
    }
    
    if (inverse) {
        uint64_t inv_n = m61_mod_pow(n, M61_MOD - 2);
        for (size_t i = 0; i < n; i++) {
            a[i] = m61_mul(a[i], inv_n);
        }
    }
    
    free(roots);
}

BigUInt *biguint_mul_ntt_mont(const BigUInt *a, const BigUInt *b) {
    if (biguint_is_zero(a) || biguint_is_zero(b)) {
        return biguint_new();
    }
    
    const uint64_t MOD = M61_MOD;
    const uint64_t ROOT = 3;
    
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
    
    ntt_inplace_m61(fa, n, ROOT, false);
    ntt_inplace_m61(fb, n, ROOT, false);
    
    for (size_t i = 0; i < n; i++) {
        fa[i] = m61_mul(fa[i], fb[i]);
    }
    
    ntt_inplace_m61(fa, n, ROOT, true);
    
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

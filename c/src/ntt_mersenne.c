#include <bigint.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#define M61_MOD (2305843009213693951ULL)
#define M61_P (61)
#define NTT_MOD (998244353ULL)
#define NTT_ROOT (3ULL)

static inline uint64_t m61_reduce(uint64_t x) {
    uint64_t t = (x >> M61_P) + (x & M61_MOD);
    t += (t >> M61_P);
    return t & M61_MOD;
}

static uint64_t mod_pow(uint64_t a, uint64_t e, uint64_t mod) {
    uint64_t res = 1;
    while (e) {
        if (e & 1) res = ( (__uint128_t)res * a ) % mod;
        a = ( (__uint128_t)a * a ) % mod;
        e >>= 1;
    }
    return res;
}

static inline uint64_t mod_mul(uint64_t a, uint64_t b, uint64_t mod) {
    return ( (__uint128_t)a * b ) % mod;
}

static void ntt_inplace(uint64_t *a, size_t n, uint64_t mod, uint64_t root, bool inverse) {
    for (size_t i = 1, j = 0; i < n; i++) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            uint64_t t = a[i]; a[i] = a[j]; a[j] = t;
        }
    }
    
    size_t len = 2;
    while (len <= n) {
        uint64_t wlen = mod_pow(root, (mod - 1) / len, mod);
        if (inverse) wlen = mod_pow(wlen, mod - 2, mod);
        
        uint64_t w = 1;
        for (size_t j = 0; j < len / 2; j++) {
            for (size_t i = 0; i < n; i += len) {
                uint64_t u = a[i + j];
                uint64_t v = mod_mul(a[i + j + len / 2], w, mod);
                a[i + j] = u + v;
                if (a[i + j] >= mod) a[i + j] -= mod;
                a[i + j + len / 2] = u + mod - v;
                if (a[i + j + len / 2] >= mod) a[i + j + len / 2] -= mod;
            }
            w = mod_mul(w, wlen, mod);
        }
        len <<= 1;
    }
    
    if (inverse) {
        uint64_t inv_n = mod_pow(n, mod - 2, mod);
        for (size_t i = 0; i < n; i++) {
            a[i] = mod_mul(a[i], inv_n, mod);
        }
    }
}

static void *xrealloc(void *ptr, size_t size) {
    void *p = realloc(ptr, size);
    if (!p && size > 0) abort();
    return p;
}

BigUInt *biguint_mul_ntt_mont(const BigUInt *a, const BigUInt *b) {
    if (biguint_is_zero(a) || biguint_is_zero(b)) {
        return biguint_new();
    }
    
    size_t n = 1;
    while (n < a->len + b->len - 1) n <<= 1;
    size_t result_len = a->len + b->len - 1;
    
    uint64_t *fa = calloc(n, sizeof(uint64_t));
    uint64_t *fb = calloc(n, sizeof(uint64_t));
    
    if (!fa || !fb) {
        free(fa); free(fb);
        return NULL;
    }
    
    for (size_t i = 0; i < a->len; i++) fa[i] = a->words[i] % NTT_MOD;
    for (size_t i = 0; i < b->len; i++) fb[i] = b->words[i] % NTT_MOD;
    
    ntt_inplace(fa, n, NTT_MOD, NTT_ROOT, false);
    ntt_inplace(fb, n, NTT_MOD, NTT_ROOT, false);
    
    for (size_t i = 0; i < n; i++) {
        fa[i] = mod_mul(fa[i], fb[i], NTT_MOD);
    }
    
    ntt_inplace(fa, n, NTT_MOD, NTT_ROOT, true);
    
    BigUInt *res = biguint_new();
    if (!res) {
        free(fa); free(fb);
        return NULL;
    }
    
    res->words = xrealloc(NULL, result_len * sizeof(uint64_t));
    if (!res->words) {
        free(fa); free(fb);
        free(res);
        return NULL;
    }
    
    for (size_t i = 0; i < result_len; i++) {
        res->words[i] = m61_reduce(fa[i]);
    }
    res->len = result_len;
    res->cap = result_len;
    
    free(fa); free(fb);
    biguint_normalize(res);
    
    return res;
}

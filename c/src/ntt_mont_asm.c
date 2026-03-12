#include <bigint.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#define NTT_MOD (998244353ULL)
#define NTT_ROOT (3ULL)
#define MONTGOMERY_R (4294967296ULL)
#define MONTGOMERY_R_INV (330382944ULL)

static inline uint64_t montgomery_mul(uint64_t a, uint64_t b, uint64_t mod) {
    uint64_t result;
    __asm__ (
        "mulq %2\n"
        "movq %%rax, %%r11\n"
        "movq %%rdx, %%r12\n"
        "mulq %3\n"
        "addq %%r11, %%rax\n"
        "adcq %%r12, %%rdx\n"
        "movq $0, %%r11\n"
        "cmpq %4, %%rdx\n"
        "cmovaeq %%r11, %%rdx\n"
        : "=a"(result)
        : "a"(a), "r"(b), "r"(MONTGOMERY_R_INV), "r"(mod)
        : "r11", "r12", "rdx", "cc"
    );
    return result;
}

static uint64_t mod_pow_asm(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t result = 1;
    uint64_t b = base % mod;
    
    while (exp > 0) {
        if (exp & 1) {
            result = montgomery_mul(result, b, mod);
        }
        b = montgomery_mul(b, b, mod);
        exp >>= 1;
    }
    return result;
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
        uint64_t wlen = mod_pow_asm(root, (mod - 1) / len, mod);
        if (inverse) wlen = mod_pow_asm(wlen, mod - 2, mod);
        
        uint64_t w = 1;
        for (size_t j = 0; j < len / 2; j++) {
            for (size_t i = 0; i < n; i += len) {
                uint64_t u = a[i + j];
                uint64_t v = montgomery_mul(a[i + j + len / 2], w, mod);
                a[i + j] = u + v;
                if (a[i + j] >= mod) a[i + j] -= mod;
                a[i + j + len / 2] = u + mod - v;
                if (a[i + j + len / 2] >= mod) a[i + j + len / 2] -= mod;
            }
            w = montgomery_mul(w, wlen, mod);
        }
        len <<= 1;
    }
    
    if (inverse) {
        uint64_t inv_n = mod_pow_asm(n, mod - 2, mod);
        for (size_t i = 0; i < n; i++) {
            a[i] = montgomery_mul(a[i], inv_n, mod);
        }
    }
}

static void *xrealloc(void *ptr, size_t size) {
    void *p = realloc(ptr, size);
    if (!p && size > 0) abort();
    return p;
}

BigUInt *biguint_mul_ntt_mont_asm(const BigUInt *a, const BigUInt *b) {
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
        fa[i] = montgomery_mul(fa[i], fb[i], NTT_MOD);
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
        res->words[i] = fa[i];
    }
    res->len = result_len;
    res->cap = result_len;
    
    free(fa); free(fb);
    biguint_normalize(res);
    
    return res;
}

#include <bigint.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#define NTT_MOD (998244353ULL)
#define NTT_ROOT (3ULL)
/* MONTGOMERY_R_INV = (-NTT_MOD)^-1 mod 2^64 */
#define MONTGOMERY_R_INV (17450252288407896063ULL)

/*
 * montgomery_mul(a, b, mod) computes (a * b * R^-1) mod `mod`
 * where R = 2^64.
 *
 * It requires `a` and `b` to be in Montgomery form: a = A*R mod mod.
 * The constant MONTGOMERY_R_INV must be (-mod)^-1 mod 2^64.
 */
static inline uint64_t montgomery_mul(uint64_t a, uint64_t b, uint64_t mod) {
    uint64_t result;
    /*
     * 1. Compute T = a * b (128-bit product in rdx:rax)
     * 2. Compute m = (T_lo * MONTGOMERY_R_INV) mod 2^64
     * 3. Compute T + m * mod. The low 64 bits are guaranteed 0.
     * 4. The result is the high 64 bits, minus mod if it overflowed >= mod.
     */
    __asm__ (
        "mulq %2\n\t"               // rdx:rax = a * b
        "movq %%rax, %%r8\n\t"        // r8 = T_lo
        "movq %%rdx, %%r9\n\t"        // r9 = T_hi
        
        "imulq %3, %%rax\n\t"         // rax = m = T_lo * mod_inv mod 2^64
        "mulq %4\n\t"                 // rdx:rax = m * mod
        
        "addq %%r8, %%rax\n\t"        // rax = T_lo + (m * mod)_lo  (always 0, sets carry)
        "adcq %%r9, %%rdx\n\t"        // rdx = T_hi + (m * mod)_hi + carry  (this is (T + m*mod)/2^64)
        
        "movq %%rdx, %%rax\n\t"       // rax = result
        "subq %4, %%rdx\n\t"          // rdx = result - mod
        "cmovnc %%rdx, %%rax\n\t"     // if (result >= mod) result = result - mod
        : "=&a"(result)               // %0
        : "0"(a),                     // %1 (same as %0, rax)
          "r"(b),                     // %2
          "r"(MONTGOMERY_R_INV),      // %3
          "r"(mod)                    // %4
        : "rdx", "r8", "r9", "cc"
    );
    return result;
}

static uint64_t mod_pow_asm(uint64_t base, uint64_t exp, uint64_t mod) {
    // 1 in Montgomery form is R mod MOD. R = 2^64. 2^64 % 998244353 = 932051910
    uint64_t result = 932051910ULL;
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

static void ntt_inplace_asm(uint64_t *a, size_t n, uint64_t mod, uint64_t root, bool inverse) {
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
        
        uint64_t w = 932051910ULL; // 1 in Montgomery form
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
        uint64_t mont_n = montgomery_mul(n % mod, 299560064ULL, mod);
        uint64_t inv_n = mod_pow_asm(mont_n, mod - 2, mod);
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
    
    /*
     * R = 2^64 mod NTT_MOD. We need this to encode inputs into Montgomery form.
     * 2^64 = 18446744073709551616.
     * 18446744073709551616 % 998244353 = 136407672
     */
    /* R^2 mod NTT_MOD. Used to convert normal -> Montgomery */
    const uint64_t R2_MOD = 299560064ULL;
    
    // Encode inputs to Montgomery form: a_mont = a * R mod NTT_MOD
    // We achieve this by montgomery_mul(a, R2_MOD) = a * R^2 * R^-1 = a * R
    for (size_t i = 0; i < a->len; i++) fa[i] = montgomery_mul(a->words[i] % NTT_MOD, R2_MOD, NTT_MOD);
    for (size_t i = 0; i < b->len; i++) fb[i] = montgomery_mul(b->words[i] % NTT_MOD, R2_MOD, NTT_MOD);
    
    // The NTT_ROOT also needs to be in Montgomery form
    uint64_t mont_root = montgomery_mul(NTT_ROOT, R2_MOD, NTT_MOD);
    
    ntt_inplace_asm(fa, n, NTT_MOD, mont_root, false);
    ntt_inplace_asm(fb, n, NTT_MOD, mont_root, false);
    
    for (size_t i = 0; i < n; i++) {
        fa[i] = montgomery_mul(fa[i], fb[i], NTT_MOD);
    }
    
    ntt_inplace_asm(fa, n, NTT_MOD, mont_root, true);
    
    // Decode from Montgomery form: montgomery_mul(x, 1) = x * 1 * R^-1 = x * R^-1
    // which cancels the R factor.
    for (size_t i = 0; i < n; i++) {
        fa[i] = montgomery_mul(fa[i], 1, NTT_MOD);
    }
    
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

#include <bigint.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#include <immintrin.h>

#define M61_MOD (2305843009213693951ULL)
#define M61_P (61)

#define M61_MASK _mm256_set1_epi64x(M61_MOD)

static inline __m256i m256i_loadu(const uint64_t *p) {
    return _mm256_loadu_si256((const __m256i *)p);
}

static inline void m256i_storeu(uint64_t *p, __m256i v) {
    _mm256_storeu_si256((__m256i *)p, v);
}

static inline __m256i m61_reduce_avx(__m256i x) {
    __m256i mod = M61_MASK;
    __m256i t = _mm256_srli_epi64(x, M61_P);
    t = _mm256_add_epi64(t, _mm256_and_si256(x, mod));
    t = _mm256_add_epi64(t, _mm256_srli_epi64(t, M61_P));
    return _mm256_and_si256(t, mod);
}

static inline __m256i m61_add_avx(__m256i a, __m256i b) {
    __m256i x = _mm256_add_epi64(a, b);
    x = _mm256_add_epi64(x, _mm256_srli_epi64(x, M61_P));
    return _mm256_and_si256(x, M61_MASK);
}

static inline __m256i m61_sub_avx(__m256i a, __m256i b) {
    const __m256i mod = M61_MASK;
    __m256i x = _mm256_add_epi64(a, mod);
    x = _mm256_sub_epi64(x, b);
    x = _mm256_add_epi64(x, _mm256_srli_epi64(x, M61_P));
    return _mm256_and_si256(x, mod);
}

static inline __m256i m61_mul_avx(__m256i a, __m256i b) {
    __m256i a_lo = _mm256_and_si256(a, M61_MASK);
    __m256i a_hi = _mm256_srli_epi64(a, M61_P);
    __m256i b_lo = _mm256_and_si256(b, M61_MASK);
    __m256i b_hi = _mm256_srli_epi64(b, M61_P);
    
    __m256i lo_lo = _mm256_mul_epu32(a_lo, b_lo);
    __m256i lo_hi = _mm256_mul_epu32(a_lo, b_hi);
    __m256i hi_lo = _mm256_mul_epu32(a_hi, b_lo);
    __m256i hi_hi = _mm256_mul_epu32(a_hi, b_hi);
    
    __m256i t1 = _mm256_or_si256(lo_lo, _mm256_slli_epi64(lo_hi, 32));
    __m256i t2 = _mm256_or_si256(hi_lo, _mm256_slli_epi64(hi_hi, 32));
    __m256i prod = _mm256_add_epi64(t1, _mm256_slli_epi64(t2, M61_P));
    
    return m61_reduce_avx(prod);
}

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

static void ntt_inplace_m61_avx(uint64_t *a, size_t n, uint64_t root, bool inverse) {
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

static void ntt_inplace_m61_avx2(uint64_t *a, size_t n, uint64_t root, bool inverse) {
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
        size_t avx_len = len / 2;
        for (size_t i = 0; i < n; i += len) {
            size_t j = 0;
            for (; j + 8 <= avx_len; j += 8) {
                __m256i u0 = m256i_loadu(&a[i + j]);
                __m256i u1 = m256i_loadu(&a[i + j + 4]);
                __m256i v0 = m61_mul_avx(m256i_loadu(&a[i + j + avx_len]), m256i_loadu(&roots[j]));
                __m256i v1 = m61_mul_avx(m256i_loadu(&a[i + j + avx_len + 4]), m256i_loadu(&roots[j + 4]));
                
                __m256i r0 = m61_add_avx(u0, v0);
                __m256i r1 = m61_add_avx(u1, v1);
                __m256i s0 = m61_sub_avx(u0, v0);
                __m256i s1 = m61_sub_avx(u1, v1);
                
                m256i_storeu(&a[i + j], r0);
                m256i_storeu(&a[i + j + 4], r1);
                m256i_storeu(&a[i + j + avx_len], s0);
                m256i_storeu(&a[i + j + avx_len + 4], s1);
            }
            for (; j < avx_len; j++) {
                uint64_t u = a[i + j];
                uint64_t v = m61_mul(a[i + j + avx_len], roots[j]);
                a[i + j] = m61_add(u, v);
                a[i + j + avx_len] = m61_sub(u, v);
            }
        }
        len <<= 1;
    }
    
    if (inverse) {
        uint64_t inv_n = m61_mod_pow(n, M61_MOD - 2);
        __m256i inv_n_vec = _mm256_set1_epi64x(inv_n);
        for (size_t i = 0; i + 4 <= n; i += 4) {
            __m256i val = m256i_loadu(&a[i]);
            val = m61_mul_avx(val, inv_n_vec);
            m256i_storeu(&a[i], val);
        }
        for (size_t i = (n & ~3ULL); i < n; i++) {
            a[i] = m61_mul(a[i], inv_n);
        }
    }
    
    free(roots);
}

BigUInt *biguint_mul_ntt_mont_avx(const BigUInt *a, const BigUInt *b) {
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
    
    ntt_inplace_m61_avx2(fa, n, ROOT, false);
    ntt_inplace_m61_avx2(fb, n, ROOT, false);
    
    for (size_t i = 0; i + 4 <= n; i += 4) {
        __m256i a_vec = m256i_loadu(&fa[i]);
        __m256i b_vec = m256i_loadu(&fb[i]);
        __m256i prod = m61_mul_avx(a_vec, b_vec);
        m256i_storeu(&fa[i], prod);
    }
    for (size_t i = (n & ~3ULL); i < n; i++) {
        fa[i] = m61_mul(fa[i], fb[i]);
    }
    
    ntt_inplace_m61_avx2(fa, n, ROOT, true);
    
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

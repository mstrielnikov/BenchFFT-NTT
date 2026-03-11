#include <bigint.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>

static void multiply_schoolbook(const uint64_t *a, size_t a_len, 
                              const uint64_t *b, size_t b_len,
                              uint64_t *result, size_t result_len) {
    memset(result, 0, result_len * sizeof(uint64_t));
    
    for (size_t i = 0; i < a_len; i++) {
        __uint128_t carry = 0;
        for (size_t j = 0; j < b_len; j++) {
            __uint128_t prod = (__uint128_t)a[i] * b[j] + result[i + j] + carry;
            result[i + j] = (uint64_t)prod;
            carry = prod >> 64;
        }
        if (carry) {
            result[i + b_len] += (uint64_t)carry;
        }
    }
}

BigUInt *biguint_mul_schoolbook(const BigUInt *a, const BigUInt *b) {
    if (biguint_is_zero(a) || biguint_is_zero(b)) {
        return biguint_new();
    }
    
    size_t result_len = a->len + b->len;
    
    BigUInt *res = biguint_new();
    if (!res) {
        return NULL;
    }
    
    res->words = calloc(result_len, sizeof(uint64_t));
    if (!res->words) {
        free(res);
        return NULL;
    }
    res->len = result_len;
    res->cap = result_len;
    
    multiply_schoolbook(a->words, a->len, b->words, b->len, res->words, result_len);
    
    biguint_normalize(res);
    return res;
}

#if defined(__AVX2__)
#include <immintrin.h>

static void multiply_schoolbook_avx2(const uint64_t *a, size_t a_len, 
                                   const uint64_t *b, size_t b_len,
                                   uint64_t *result, size_t result_len) {
    memset(result, 0, result_len * sizeof(uint64_t));
    
    for (size_t i = 0; i < a_len; i++) {
        __uint128_t carry = 0;
        
        size_t j = 0;
        for (; j + 4 <= b_len; j += 4) {
            __m256i b_vec = _mm256_loadu_si256((__m256i*)&b[j]);
            __m256i r_vec = _mm256_loadu_si256((__m256i*)&result[i + j]);
            
            __m256i a_vec = _mm256_set1_epi64x(a[i]);
            
            __m256i prod = _mm256_mul_epu32(a_vec, b_vec);
            
            __m256i sum = _mm256_add_epi64(prod, r_vec);
            
            _mm256_storeu_si256((__m256i*)&result[i + j], sum);
            
            carry = 0;
        }
        
        for (; j < b_len; j++) {
            __uint128_t prod = (__uint128_t)a[i] * b[j] + result[i + j] + carry;
            result[i + j] = (uint64_t)prod;
            carry = prod >> 64;
        }
        
        if (carry) {
            result[i + b_len] += (uint64_t)carry;
        }
    }
}

BigUInt *biguint_mul_schoolbook_avx(const BigUInt *a, const BigUInt *b) {
    if (biguint_is_zero(a) || biguint_is_zero(b)) {
        return biguint_new();
    }
    
    size_t result_len = a->len + b->len;
    
    BigUInt *res = biguint_new();
    if (!res) {
        return NULL;
    }
    
    res->words = calloc(result_len, sizeof(uint64_t));
    if (!res->words) {
        free(res);
        return NULL;
    }
    res->len = result_len;
    res->cap = result_len;
    
    multiply_schoolbook_avx2(a->words, a->len, b->words, b->len, res->words, result_len);
    
    biguint_normalize(res);
    return res;
}
#endif

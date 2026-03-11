#include <bigint.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <immintrin.h>

static void multiply_base256_avx(const uint64_t *a, size_t a_len, 
                                 const uint64_t *b, size_t b_len,
                                 uint64_t *result, size_t result_len) {
    memset(result, 0, result_len * sizeof(uint64_t));
    
    for (size_t i = 0; i < a_len; i++) {
        uint64_t ai = a[i];
        
        for (size_t j = 0; j < b_len; j++) {
            uint64_t bi = b[j];
            
            uint8_t *a_bytes = (uint8_t *)&ai;
            uint8_t *b_bytes = (uint8_t *)&bi;
            
            for (size_t k = 0; k < 8; k++) {
                for (size_t l = 0; l < 8; l++) {
                    size_t idx = i + j + k + l;
                    if (idx < result_len) {
                        result[idx] += (uint64_t)a_bytes[k] * (uint64_t)b_bytes[l];
                    }
                }
            }
        }
    }
    
    uint64_t carry = 0;
    for (size_t i = 0; i < result_len; i++) {
        uint64_t val = result[i] + carry;
        result[i] = val & 0xFFFFFFFFFFFFFFFFULL;
        carry = val >> 64;
    }
    
    size_t block = 0;
    while (block < result_len) {
        uint64_t block_sum = 0;
        size_t block_end = block + 4;
        if (block_end > result_len) block_end = result_len;
        
        for (size_t i = block; i < block_end; i++) {
            uint64_t part = result[i] >> 32;
            result[i] = (result[i] & 0xFFFFFFFFULL) | (part << 32);
            block_sum += part;
        }
        
        uint64_t carry32 = block_sum >> 32;
        for (size_t i = block; i < block_end; i++) {
            result[i] += carry32;
            carry32 = result[i] >> 64;
            result[i] &= 0xFFFFFFFFFFFFFFFFULL;
        }
        
        block = block_end;
    }
    
    for (size_t i = 0; i < result_len - 1; i++) {
        uint64_t carry = result[i] >> 32;
        result[i] &= 0xFFFFFFFFULL;
        result[i + 1] += carry;
    }
    
    uint64_t carry = result[result_len - 1] >> 32;
    result[result_len - 1] &= 0xFFFFFFFFULL;
    
    for (size_t i = 0; i < result_len && carry > 0; i++) {
        uint64_t sum = result[i] + (carry & 0xFFFFFFFFULL);
        result[i] = sum & 0xFFFFFFFFFFFFFFFFULL;
        carry = sum >> 64;
        carry += (carry >> 32);
    }
}

BigUInt *biguint_mul_base256(const BigUInt *a, const BigUInt *b) {
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
    
    multiply_base256_avx(a->words, a->len, b->words, b->len, res->words, result_len);
    
    biguint_normalize(res);
    return res;
}

static void multiply_base256_scalar(const uint64_t *a, size_t a_len, 
                                    const uint64_t *b, size_t b_len,
                                    uint64_t *result, size_t result_len) {
    memset(result, 0, result_len * sizeof(uint64_t));
    
    for (size_t i = 0; i < a_len; i++) {
        for (size_t j = 0; j < b_len; j++) {
            __uint128_t prod = (__uint128_t)a[i] * b[j];
            result[i + j] += (uint64_t)prod;
            result[i + j + 1] += (uint64_t)(prod >> 64);
        }
    }
    
    for (size_t i = 1; i < result_len; i++) {
        result[i] += result[i - 1] >> 64;
        result[i - 1] &= 0xFFFFFFFFFFFFFFFFULL;
    }
}

BigUInt *biguint_mul_base256_scalar(const BigUInt *a, const BigUInt *b) {
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
    
    multiply_base256_scalar(a->words, a->len, b->words, b->len, res->words, result_len);
    
    biguint_normalize(res);
    return res;
}

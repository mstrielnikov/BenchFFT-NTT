#include <bigint.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>


static void multiply_standard(const uint64_t *a, size_t a_len, 
                              const uint64_t *b, size_t b_len,
                              uint64_t *result, size_t result_len) {
    memset(result, 0, result_len * sizeof(uint64_t));
    
    for (size_t i = 0; i < a_len; i++) {
        __uint128_t carry = 0;
        for (size_t j = 0; j < b_len; j++) {
            if (i + j < result_len) {
                __uint128_t prod = (__uint128_t)a[i] * b[j] + result[i + j] + carry;
                result[i + j] = (uint64_t)prod;
                carry = prod >> 64;
            }
        }
        if (carry && i + b_len < result_len) {
            result[i + b_len] += (uint64_t)carry;
        }
    }
}

BigUInt *biguint_mul_standard(const BigUInt *a, const BigUInt *b) {
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
    
    multiply_standard(a->words, a->len, b->words, b->len, res->words, result_len);
    
    biguint_normalize(res);
    return res;
}

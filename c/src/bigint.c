#include <bigint.h>
#include <stdlib.h>
#include <string.h>

static void *xrealloc(void *ptr, size_t size) {
    void *p = realloc(ptr, size);
    if (!p && size > 0) abort();
    return p;
}

void biguint_normalize(BigUInt *a) {
    while (a->len > 0 && a->words[a->len - 1] == 0) {
        a->len--;
    }
}

BigUInt *biguint_new(void) {
    BigUInt *a = malloc(sizeof(BigUInt));
    if (!a) return NULL;
    a->words = NULL;
    a->len = 0;
    a->cap = 0;
    return a;
}

BigUInt *biguint_from_uint64(uint64_t n) {
    BigUInt *a = biguint_new();
    if (!a) return NULL;
    if (n != 0) {
        a->words = malloc(sizeof(uint64_t));
        if (!a->words) { free(a); return NULL; }
        a->words[0] = n;
        a->len = 1;
        a->cap = 1;
    }
    return a;
}

BigUInt *biguint_from_slice(const uint64_t *words, size_t len) {
    BigUInt *a = biguint_new();
    if (!a) return NULL;
    if (len == 0) return a;
    
    a->words = malloc(len * sizeof(uint64_t));
    if (!a->words) { free(a); return NULL; }
    
    memcpy(a->words, words, len * sizeof(uint64_t));
    a->len = len;
    a->cap = len;
    biguint_normalize(a);
    return a;
}

BigUInt *biguint_clone(const BigUInt *a) {
    return biguint_from_slice(a->words, a->len);
}

void biguint_free(BigUInt *a) {
    if (a) {
        free(a->words);
        free(a);
    }
}

bool biguint_is_zero(const BigUInt *a) {
    return a->len == 0;
}

int biguint_cmp(const BigUInt *a, const BigUInt *b) {
    if (a->len != b->len) {
        return a->len < b->len ? -1 : 1;
    }
    for (size_t i = a->len; i > 0; i--) {
        if (a->words[i-1] != b->words[i-1]) {
            return a->words[i-1] < b->words[i-1] ? -1 : 1;
        }
    }
    return 0;
}

BigUInt *biguint_add(const BigUInt *a, const BigUInt *b) {
    size_t max_len = a->len > b->len ? a->len : b->len;
    BigUInt *res = biguint_new();
    if (!res) return NULL;
    
    res->words = calloc(max_len + 1, sizeof(uint64_t));
    if (!res->words) { biguint_free(res); return NULL; }
    res->cap = max_len + 1;
    
    for (size_t i = 0; i < max_len; i++) {
        uint64_t av = i < a->len ? a->words[i] : 0;
        uint64_t bv = i < b->len ? b->words[i] : 0;
        __uint128_t sum = (__uint128_t)res->words[i] + av + bv;
        res->words[i] = (uint64_t)sum;
        res->words[i + 1] += (uint64_t)(sum >> 64);
    }
    
    res->len = max_len + 1;
    biguint_normalize(res);
    return res;
}

size_t next_power_of_two(size_t n) {
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n |= n >> 32;
    return n + 1;
}

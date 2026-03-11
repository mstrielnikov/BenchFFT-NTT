#ifndef BIGINT_H
#define BIGINT_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#if defined(__AVX2__) || defined(__AVX__)
#define HAS_AVX 1
#else
#define HAS_AVX 0
#endif

#if defined(BUILD_AVX) && HAS_AVX
#define biguint_mul_fft_split biguint_mul_fft_split_avx
#define biguint_mul_ntt_mont biguint_mul_ntt_mont_avx
#endif

typedef struct {
    uint64_t *words;
    size_t len;
    size_t cap;
} BigUInt;

BigUInt *biguint_new(void);
BigUInt *biguint_from_uint64(uint64_t n);
BigUInt *biguint_from_slice(const uint64_t *words, size_t len);
BigUInt *biguint_clone(const BigUInt *a);
void biguint_free(BigUInt *a);

bool biguint_is_zero(const BigUInt *a);
void biguint_normalize(BigUInt *a);
int biguint_cmp(const BigUInt *a, const BigUInt *b);

BigUInt *biguint_add(const BigUInt *a, const BigUInt *b);

BigUInt *biguint_mul_fft_split(const BigUInt *a, const BigUInt *b);
BigUInt *biguint_mul_ntt_mont(const BigUInt *a, const BigUInt *b);

#if HAS_AVX
BigUInt *biguint_mul_fft_split_avx(const BigUInt *a, const BigUInt *b);
BigUInt *biguint_mul_ntt_mont_avx(const BigUInt *a, const BigUInt *b);
#endif

size_t next_power_of_two(size_t n);

#endif

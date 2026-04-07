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

typedef struct {
    uint64_t *words;
    size_t len;
    size_t cap;
} BigUInt;

void biguint_free(BigUInt *a);
bool biguint_is_zero(const BigUInt *a);
void biguint_normalize(BigUInt *a);
int biguint_cmp(const BigUInt *a, const BigUInt *b);
size_t biguint_len(const BigUInt *a);
uint64_t biguint_get_word(const BigUInt *a, size_t i);
size_t next_power_of_two(size_t n);

BigUInt *biguint_new(void);
BigUInt *biguint_from_uint64(uint64_t n);
BigUInt *biguint_from_slice(const uint64_t *words, size_t len);

BigUInt *biguint_mul_standard(const BigUInt *a, const BigUInt *b);

BigUInt *biguint_mul_fft_split(const BigUInt *a, const BigUInt *b);
BigUInt *biguint_mul_fft_mersenne(const BigUInt *a, const BigUInt *b);
BigUInt *biguint_mul_ntt_mersenne(const BigUInt *a, const BigUInt *b);
BigUInt *biguint_mul_ntt_mont(const BigUInt *a, const BigUInt *b);
BigUInt *biguint_mul_ntt_mont_m61(const BigUInt *a, const BigUInt *b);   /* Integer-NTT mod 998244353 → M61 reduction*/
BigUInt *biguint_mul_ntt_crt(const BigUInt *a, const BigUInt *b);        /* Multi-modular NTT + CRT exact reconstruction */
BigUInt *biguint_mul_nussbaumer(const BigUInt *a, const BigUInt *b);
BigUInt *biguint_mul_bluestein(const BigUInt *a, const BigUInt *b);
BigUInt *biguint_mul_toom3(const BigUInt *a, const BigUInt *b);


#if HAS_AVX
BigUInt *biguint_mul_standard_avx(const BigUInt *a, const BigUInt *b);
BigUInt *biguint_mul_fft_split_avx(const BigUInt *a, const BigUInt *b);
BigUInt *biguint_mul_ntt_mersenne_avx(const BigUInt *a, const BigUInt *b);
BigUInt *biguint_mul_ntt_mont_avx(const BigUInt *a, const BigUInt *b);
BigUInt *biguint_mul_ntt_mont_m61_avx(const BigUInt *a, const BigUInt *b);
BigUInt *biguint_mul_toom3_avx(const BigUInt *a, const BigUInt *b);
#endif

#ifdef __x86_64__
BigUInt *biguint_mul_ntt_mont_asm(const BigUInt *a, const BigUInt *b);
#endif

#endif

#include <bigint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int tests_passed = 0;
static int tests_failed = 0;

#define ASSERT(cond, msg) do { \
    if (cond) { \
        tests_passed++; \
    } else { \
        tests_failed++; \
        printf("FAIL: %s\n", msg); \
    } \
} while(0)

#define ASSERT_EQ(a, b, msg) do { \
    int cmp = biguint_cmp(a, b); \
    if (cmp == 0) { \
        tests_passed++; \
    } else { \
        tests_failed++; \
        printf("FAIL: %s\n", msg); \
    } \
} while(0)

void test_fft_split() {
    printf("=== FFT SPLIT TESTS ===\n");
    
    // Test small multiplication
    BigUInt *a = biguint_from_uint64(12345);
    BigUInt *b = biguint_from_uint64(67890);
    BigUInt *c = biguint_mul_fft_split(a, b);
    BigUInt *expected = biguint_from_uint64(12345ULL * 67890ULL);
    ASSERT_EQ(c, expected, "12345 * 67890");
    biguint_free(c);
    biguint_free(expected);
    biguint_free(a);
    biguint_free(b);
    
    // Test ML-KEM vector size: 256 words
    uint64_t *words256_1 = malloc(256 * sizeof(uint64_t));
    uint64_t *words256_2 = malloc(256 * sizeof(uint64_t));
    for (int i = 0; i < 256; i++) {
        words256_1[i] = i + 1;
        words256_2[i] = (i + 1) * 2;
    }
    a = biguint_from_slice(words256_1, 256);
    b = biguint_from_slice(words256_2, 256);
    c = biguint_mul_fft_split(a, b);
    ASSERT(c != NULL && c->len > 256, "mul 256 words");
    biguint_free(c);
    biguint_free(a);
    biguint_free(b);
    free(words256_1);
    free(words256_2);
    
    // Test ML-KEM vector size: 512 words
    uint64_t *words512_1 = malloc(512 * sizeof(uint64_t));
    uint64_t *words512_2 = malloc(512 * sizeof(uint64_t));
    for (int i = 0; i < 512; i++) {
        words512_1[i] = i + 1;
        words512_2[i] = (i + 1) * 2;
    }
    a = biguint_from_slice(words512_1, 512);
    b = biguint_from_slice(words512_2, 512);
    c = biguint_mul_fft_split(a, b);
    ASSERT(c != NULL && c->len > 512, "mul 512 words");
    biguint_free(c);
    biguint_free(a);
    biguint_free(b);
    free(words512_1);
    free(words512_2);
    
    // Test ML-KEM vector size: 1024 words
    uint64_t *words1024_1 = malloc(1024 * sizeof(uint64_t));
    uint64_t *words1024_2 = malloc(1024 * sizeof(uint64_t));
    for (int i = 0; i < 1024; i++) {
        words1024_1[i] = i + 1;
        words1024_2[i] = (i + 1) * 2;
    }
    a = biguint_from_slice(words1024_1, 1024);
    b = biguint_from_slice(words1024_2, 1024);
    c = biguint_mul_fft_split(a, b);
    ASSERT(c != NULL && c->len > 1024, "mul 1024 words");
    biguint_free(c);
    biguint_free(a);
    biguint_free(b);
    free(words1024_1);
    free(words1024_2);
    
    // Test ML-DSA vector size: 2048 words
    uint64_t *words2048_1 = malloc(2048 * sizeof(uint64_t));
    uint64_t *words2048_2 = malloc(2048 * sizeof(uint64_t));
    for (int i = 0; i < 2048; i++) {
        words2048_1[i] = i + 1;
        words2048_2[i] = (i + 1) * 2;
    }
    a = biguint_from_slice(words2048_1, 2048);
    b = biguint_from_slice(words2048_2, 2048);
    c = biguint_mul_fft_split(a, b);
    ASSERT(c != NULL && c->len > 2048, "mul 2048 words");
    biguint_free(c);
    biguint_free(a);
    biguint_free(b);
    free(words2048_1);
    free(words2048_2);
    
    // Test Extended: 3072 words
    uint64_t *words3072_1 = malloc(3072 * sizeof(uint64_t));
    uint64_t *words3072_2 = malloc(3072 * sizeof(uint64_t));
    for (int i = 0; i < 3072; i++) {
        words3072_1[i] = i + 1;
        words3072_2[i] = (i + 1) * 2;
    }
    a = biguint_from_slice(words3072_1, 3072);
    b = biguint_from_slice(words3072_2, 3072);
    c = biguint_mul_fft_split(a, b);
    ASSERT(c != NULL && c->len > 3072, "mul 3072 words");
    biguint_free(c);
    biguint_free(a);
    biguint_free(b);
    free(words3072_1);
    free(words3072_2);
    
    // Test Extended: 4096 words
    uint64_t *words4096_1 = malloc(4096 * sizeof(uint64_t));
    uint64_t *words4096_2 = malloc(4096 * sizeof(uint64_t));
    for (int i = 0; i < 4096; i++) {
        words4096_1[i] = i + 1;
        words4096_2[i] = (i + 1) * 2;
    }
    a = biguint_from_slice(words4096_1, 4096);
    b = biguint_from_slice(words4096_2, 4096);
    c = biguint_mul_fft_split(a, b);
    ASSERT(c != NULL && c->len > 4096, "mul 4096 words");
    biguint_free(c);
    biguint_free(a);
    biguint_free(b);
    free(words4096_1);
    free(words4096_2);
    
    // Test with max values
    uint64_t max_words[] = {0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF};
    a = biguint_from_slice(max_words, 2);
    b = biguint_from_slice(max_words, 2);
    c = biguint_mul_fft_split(a, b);
    ASSERT(c != NULL, "mul max values");
    biguint_free(c);
    biguint_free(a);
    biguint_free(b);
    
    printf("FFT Split: %d passed, %d failed\n", tests_passed, tests_failed);
}

void test_ntt_mont() {
    printf("=== NTT MONT TESTS ===\n");
    int start_passed = tests_passed;
    
    // Test small multiplication
    BigUInt *a = biguint_from_uint64(12345);
    BigUInt *b = biguint_from_uint64(67890);
    BigUInt *c = biguint_mul_ntt_mont(a, b);
    BigUInt *expected = biguint_from_uint64(12345ULL * 67890ULL);
    ASSERT_EQ(c, expected, "12345 * 67890");
    biguint_free(c);
    biguint_free(expected);
    biguint_free(a);
    biguint_free(b);
    
    // Test ML-KEM vector size: 256 words
    uint64_t *words256_1 = malloc(256 * sizeof(uint64_t));
    uint64_t *words256_2 = malloc(256 * sizeof(uint64_t));
    for (int i = 0; i < 256; i++) {
        words256_1[i] = i + 1;
        words256_2[i] = (i + 1) * 2;
    }
    a = biguint_from_slice(words256_1, 256);
    b = biguint_from_slice(words256_2, 256);
    c = biguint_mul_ntt_mont(a, b);
    ASSERT(c != NULL && c->len > 256, "mul 256 words");
    biguint_free(c);
    biguint_free(a);
    biguint_free(b);
    free(words256_1);
    free(words256_2);
    
    // Test ML-KEM vector size: 512 words
    uint64_t *words512_1 = malloc(512 * sizeof(uint64_t));
    uint64_t *words512_2 = malloc(512 * sizeof(uint64_t));
    for (int i = 0; i < 512; i++) {
        words512_1[i] = i + 1;
        words512_2[i] = (i + 1) * 2;
    }
    a = biguint_from_slice(words512_1, 512);
    b = biguint_from_slice(words512_2, 512);
    c = biguint_mul_ntt_mont(a, b);
    ASSERT(c != NULL && c->len > 512, "mul 512 words");
    biguint_free(c);
    biguint_free(a);
    biguint_free(b);
    free(words512_1);
    free(words512_2);
    
    // Test Extended: 3072 words
    uint64_t *words3072_1 = malloc(3072 * sizeof(uint64_t));
    uint64_t *words3072_2 = malloc(3072 * sizeof(uint64_t));
    for (int i = 0; i < 3072; i++) {
        words3072_1[i] = i + 1;
        words3072_2[i] = (i + 1) * 2;
    }
    a = biguint_from_slice(words3072_1, 3072);
    b = biguint_from_slice(words3072_2, 3072);
    c = biguint_mul_ntt_mont(a, b);
    ASSERT(c != NULL && c->len > 3072, "mul 3072 words");
    biguint_free(c);
    biguint_free(a);
    biguint_free(b);
    free(words3072_1);
    free(words3072_2);
    
    // Test Extended: 4096 words
    uint64_t *words4096_1 = malloc(4096 * sizeof(uint64_t));
    uint64_t *words4096_2 = malloc(4096 * sizeof(uint64_t));
    for (int i = 0; i < 4096; i++) {
        words4096_1[i] = i + 1;
        words4096_2[i] = (i + 1) * 2;
    }
    a = biguint_from_slice(words4096_1, 4096);
    b = biguint_from_slice(words4096_2, 4096);
    c = biguint_mul_ntt_mont(a, b);
    ASSERT(c != NULL && c->len > 4096, "mul 4096 words");
    biguint_free(c);
    biguint_free(a);
    biguint_free(b);
    free(words4096_1);
    free(words4096_2);
    
    printf("NTT Mont: %d passed, %d failed\n", tests_passed - start_passed, tests_failed);
}

void test_add() {
    printf("=== ADD TESTS ===\n");
    int start_passed = tests_passed;
    
    BigUInt *a = biguint_from_uint64(10);
    BigUInt *b = biguint_from_uint64(20);
    BigUInt *c = biguint_add(a, b);
    BigUInt *expected = biguint_from_uint64(30);
    ASSERT_EQ(c, expected, "10 + 20");
    biguint_free(c);
    biguint_free(expected);
    biguint_free(a);
    biguint_free(b);
    
    printf("Add: %d passed, %d failed\n", tests_passed - start_passed, tests_failed);
}

int main() {
    test_add();
    test_fft_split();
    test_ntt_mont();
    
    printf("\n====================\n");
    printf("TOTAL: %d passed, %d failed\n", tests_passed, tests_failed);
    printf("====================\n");
    return tests_failed > 0 ? 1 : 0;
}

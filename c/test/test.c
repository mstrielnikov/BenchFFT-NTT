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
    
    BigUInt *a = biguint_from_uint64(12345);
    BigUInt *b = biguint_from_uint64(67890);
    BigUInt *c = biguint_mul_fft_split(a, b);
    BigUInt *expected = biguint_from_uint64(12345ULL * 67890ULL);
    ASSERT_EQ(c, expected, "12345 * 67890");
    biguint_free(c);
    biguint_free(expected);
    biguint_free(a);
    biguint_free(b);
    
    uint64_t words1[] = {0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF};
    uint64_t words2[] = {2, 0};
    a = biguint_from_slice(words1, 2);
    b = biguint_from_slice(words2, 2);
    c = biguint_mul_fft_split(a, b);
    ASSERT(c != NULL, "mul max values");
    biguint_free(c);
    biguint_free(a);
    biguint_free(b);
    
    printf("FFT Split: %d passed, %d failed\n", tests_passed, tests_failed);
}

void test_ntt_mont() {
    printf("=== NTT MONT TESTS ===\n");
    
    BigUInt *a = biguint_from_uint64(12345);
    BigUInt *b = biguint_from_uint64(67890);
    BigUInt *c = biguint_mul_ntt_mont(a, b);
    BigUInt *expected = biguint_from_uint64(12345ULL * 67890ULL);
    ASSERT_EQ(c, expected, "12345 * 67890");
    biguint_free(c);
    biguint_free(expected);
    biguint_free(a);
    biguint_free(b);
    
    printf("NTT Mont: %d passed, %d failed\n", tests_passed - 2, tests_failed);
}

void test_add() {
    printf("=== ADD TESTS ===\n");
    
    BigUInt *a = biguint_from_uint64(10);
    BigUInt *b = biguint_from_uint64(20);
    BigUInt *c = biguint_add(a, b);
    BigUInt *expected = biguint_from_uint64(30);
    ASSERT_EQ(c, expected, "10 + 20");
    biguint_free(c);
    biguint_free(expected);
    biguint_free(a);
    biguint_free(b);
    
    printf("Add: %d passed, %d failed\n", tests_passed - 4, tests_failed);
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

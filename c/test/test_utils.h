/**
 * test_harness.h — shared test infrastructure for BenchFFT-NTT
 *
 * Include this exactly once from a single translation unit.
 * It defines the shared pass/fail counters, assertion macros, and a generic
 * run_mul_tests() that runs the standard ML-KEM/ML-DSA size battery against
 * any multiplication function pointer.
 */
#pragma once

#include <bigint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ── Shared counters ──────────────────────────────────────────────────────── */

static int tests_passed = 0;
static int tests_failed = 0;

/* ── Assertion macros ─────────────────────────────────────────────────────── */

#define ASSERT(cond, msg) do {          \
    if (cond) {                         \
        tests_passed++;                 \
    } else {                            \
        tests_failed++;                 \
        printf("FAIL: %s\n", msg);      \
    }                                   \
} while (0)

#define ASSERT_EQ(a, b, msg) do {       \
    if (biguint_cmp(a, b) == 0) {       \
        tests_passed++;                 \
    } else {                            \
        tests_failed++;                 \
        printf("FAIL: %s\n", msg);      \
    }                                   \
} while (0)

/* ── Standard test sizes ──────────────────────────────────────────────────── */

static const size_t TEST_SIZES[] = { 256, 512, 1024, 2048, 3072, 4096 };
#define N_TEST_SIZES (sizeof(TEST_SIZES) / sizeof(TEST_SIZES[0]))

/* ── Generic test battery ─────────────────────────────────────────────────── */

/**
 * run_mul_tests — run the standard multiplication correctness battery.
 *
 * @param name  Label printed in the section header
 * @param mul   Multiplication function under test
 *
 * Tests:
 *   1. Small known value: 12345 × 67890 = 838102050
 *   2. Output length > input length for each ML-KEM/ML-DSA size
 *   3. Non-NULL result for two max-value (all-ones) 64-bit words
 */
static void run_mul_tests(const char *name,
                          BigUInt *(*mul)(const BigUInt *, const BigUInt *))
{
    printf("=== %s TESTS ===\n", name);
    int pass_before = tests_passed;

    /* 1. Known small value */
    BigUInt *a = biguint_from_uint64(12345);
    BigUInt *b = biguint_from_uint64(67890);
    BigUInt *c = mul(a, b);
    BigUInt *expected = biguint_from_uint64(12345ULL * 67890ULL);
    ASSERT_EQ(c, expected, "12345 * 67890");
    biguint_free(c); biguint_free(expected);
    biguint_free(a); biguint_free(b);

    /* 2. Length battery */
    for (size_t s = 0; s < N_TEST_SIZES; s++) {
        size_t n = TEST_SIZES[s];
        uint64_t *w1 = malloc(n * sizeof(uint64_t));
        uint64_t *w2 = malloc(n * sizeof(uint64_t));
        for (size_t i = 0; i < n; i++) {
            w1[i] = i + 1;
            w2[i] = (i + 1) * 2;
        }
        a = biguint_from_slice(w1, n);
        b = biguint_from_slice(w2, n);
        c = mul(a, b);
        char msg[64];
        snprintf(msg, sizeof(msg), "mul %zu words", n);
        ASSERT(c != NULL && c->len > n, msg);
        biguint_free(c); biguint_free(a); biguint_free(b);
        free(w1); free(w2);
    }

    /* 3. Max-value sentinel */
    uint64_t max_words[] = { 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF };
    a = biguint_from_slice(max_words, 2);
    b = biguint_from_slice(max_words, 2);
    c = mul(a, b);
    ASSERT(c != NULL, "mul max values");
    biguint_free(c); biguint_free(a); biguint_free(b);

    int suite_pass = tests_passed - pass_before;
    int suite_fail = tests_failed;  /* cumulative; warn if > 0 */
    printf("%s: %d passed, %d failed\n\n", name, suite_pass, suite_fail);
}

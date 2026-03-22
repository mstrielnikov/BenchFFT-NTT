#include <bigint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <stdio.h>
#include "../src/fft_split.c"

/* ── Standard benchmark configuration ───────────────────────────────────── */

typedef struct {
    size_t size;
    int    iters;
    const char *label;
} BenchCase;

static const BenchCase BENCH_CASES[] = {
    { 256,  100, "ML-KEM-512  (256 words)" },
    { 512,   50, "ML-KEM-768  (512 words)" },
    { 1024,  20, "ML-KEM-1024 (1024 words)" },
    { 2048,  10, "ML-DSA      (2048 words)" },
    { 3072,   5, "Extended    (3072 words)" },
    { 4096,   3, "Extended    (4096 words)" },
};
#define N_BENCH_CASES (sizeof(BENCH_CASES) / sizeof(BENCH_CASES[0]))

static double time_ms(clock_t start, clock_t end) {
    return (double)(end - start) * 1000.0 / CLOCKS_PER_SEC;
}

/* ── Generic benchmark runner ────────────────────────────────────────────── */

/**
 * run_benchmark — time a multiplication function over all standard sizes.
 *
 * @param tag   Short label printed per row (e.g. "FFT", "FFT-AVX")
 * @param mul   Multiplication function under test
 */
static void run_benchmark(const char *tag,
                          BigUInt *(*mul)(const BigUInt *, const BigUInt *))
{
    for (size_t i = 0; i < N_BENCH_CASES; i++) {
        size_t n    = BENCH_CASES[i].size;
        int    iter = BENCH_CASES[i].iters;

        printf("%s:\n", BENCH_CASES[i].label);
        uint64_t *w1 = malloc(n * sizeof(uint64_t));
        uint64_t *w2 = malloc(n * sizeof(uint64_t));
        for (size_t j = 0; j < n; j++) {
            w1[j] = ((uint64_t)rand() << 32) | (uint64_t)rand();
            w2[j] = ((uint64_t)rand() << 32) | (uint64_t)rand();
        }
        BigUInt *a = biguint_from_slice(w1, n);
        BigUInt *b = biguint_from_slice(w2, n);
        free(w1); free(w2);

        clock_t start = clock();
        for (int k = 0; k < iter; k++) {
            BigUInt *c = mul(a, b);
            biguint_free(c);
        }
        clock_t end = clock();

        printf("%-9s %4zu:   %.4f ms (%d iters)\n",
               tag, n, time_ms(start, end), iter);

        biguint_free(a); biguint_free(b);
    }
}

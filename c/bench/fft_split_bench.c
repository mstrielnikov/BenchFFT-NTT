#include <bigint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include "../src/fft_split.c"

static double time_ms(clock_t start, clock_t end) {
    return (double)(end - start) * 1000.0 / CLOCKS_PER_SEC;
}

void benchmark_fft_split(const char *name, size_t size, int iterations) {
    uint64_t *words1 = malloc(size * sizeof(uint64_t));
    uint64_t *words2 = malloc(size * sizeof(uint64_t));
    for (size_t i = 0; i < size; i++) {
        words1[i] = ((uint64_t)rand() << 32) | rand();
        words2[i] = ((uint64_t)rand() << 32) | rand();
    }
    
    BigUInt *a = biguint_from_slice(words1, size);
    BigUInt *b = biguint_from_slice(words2, size);
    free(words1);
    free(words2);
    
    // FFT Benchmark
    clock_t start = clock();
    for (int i = 0; i < iterations; i++) {
        BigUInt *c = biguint_mul_fft_split(a, b);
        biguint_free(c);
    }
    clock_t end = clock();
    printf("FFT       %4zu:   %.4f ms (%d iters)\n", size, time_ms(start, end), iterations);
    
    biguint_free(a);
    biguint_free(b);
}

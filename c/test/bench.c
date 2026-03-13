#include <bigint.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

static double time_ms(clock_t start, clock_t end) {
    return (double)(end - start) * 1000.0 / CLOCKS_PER_SEC;
}

void run_benchmark(const char *name, size_t size, int iterations) {
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
    
    // FFT AVX Benchmark
    start = clock();
    for (int i = 0; i < iterations; i++) {
        BigUInt *c = biguint_mul_fft_split_avx(a, b);
        biguint_free(c);
    }
    end = clock();
    // FFT M61 Benchmark
    start = clock();
    for (int i = 0; i < iterations; i++) {
        BigUInt *c = biguint_mul_fft_mersenne(a, b);
        biguint_free(c);
    }
    end = clock();
    printf("FFT M61   %4zu:   %.4f ms (%d iters)\n", size, time_ms(start, end), iterations);
    
    // NTT M61 (Mersenne scalar) Benchmark
    start = clock();
    for (int i = 0; i < iterations; i++) {
        BigUInt *c = biguint_mul_ntt_mont_m61(a, b);
        biguint_free(c);
    }
    end = clock();
    printf("NTT M61   %4zu:   %.4f ms (%d iters)\n", size, time_ms(start, end), iterations);
    
    // NTT Benchmark
    start = clock();
    for (int i = 0; i < iterations; i++) {
        BigUInt *c = biguint_mul_ntt_mont(a, b);
        biguint_free(c);
    }
    end = clock();
    printf("NTT       %4zu:   %.4f ms (%d iters)\n", size, time_ms(start, end), iterations);
    
    // NTT ASM Benchmark
    start = clock();
    for (int i = 0; i < iterations; i++) {
        BigUInt *c = biguint_mul_ntt_mont_asm(a, b);
        biguint_free(c);
    }
    end = clock();
    printf("NTT ASM   %4zu:   %.4f ms (%d iters)\n", size, time_ms(start, end), iterations);
    
    // Mersenne (Direct M61) Benchmark
    start = clock();
    for (int i = 0; i < iterations; i++) {
        BigUInt *c = biguint_mul_ntt_mersenne_avx(a, b);
        biguint_free(c);
    }
    end = clock();
    printf("Mersenne  %4zu:   %.4f ms (%d iters)\n", size, time_ms(start, end), iterations);
    
    biguint_free(a);
    biguint_free(b);
}

int main() {
    printf("========== BENCHMARKS ==========\n\n");
    printf("ML-KEM / ML-DSA Vector Sizes\n");
    printf("=============================\n\n");
    
    srand(42);
    
    printf("ML-KEM-512 (256 words):\n");
    run_benchmark("FFT", 256, 100);
    printf("\n");
    
    printf("ML-KEM-768 (512 words):\n");
    run_benchmark("FFT", 512, 50);
    printf("\n");
    
    printf("ML-KEM-1024 (1024 words):\n");
    run_benchmark("FFT", 1024, 20);
    printf("\n");
    
    printf("ML-DSA (2048 words):\n");
    run_benchmark("FFT", 2048, 10);
    printf("\n");
    
    printf("Extended (3072 words):\n");
    run_benchmark("FFT", 3072, 5);
    printf("\n");
    
    printf("Extended (4096 words):\n");
    run_benchmark("FFT", 4096, 3);
    printf("\n");
    
    printf("================================\n");
    return 0;
}

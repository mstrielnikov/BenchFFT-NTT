#include <bigint.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

static double time_ms(clock_t start, clock_t end) {
    return (double)(end - start) * 1000.0 / CLOCKS_PER_SEC;
}

int main() {
    printf("========== BENCHMARKS ==========\n\n");
    
    srand(42);
    
    printf("=== MUL FFT variants (256 words) ===\n");
    
    uint64_t *words1 = malloc(256 * sizeof(uint64_t));
    uint64_t *words2 = malloc(256 * sizeof(uint64_t));
    for (int i = 0; i < 256; i++) {
        words1[i] = ((uint64_t)rand() << 32) | rand();
        words2[i] = ((uint64_t)rand() << 32) | rand();
    }
    
    BigUInt *a = biguint_from_slice(words1, 256);
    BigUInt *b = biguint_from_slice(words2, 256);
    free(words1);
    free(words2);
    
    clock_t start = clock();
    for (int i = 0; i < 100; i++) {
        BigUInt *c = biguint_mul_fft_split(a, b);
        biguint_free(c);
    }
    clock_t end = clock();
    printf("mul (FFT Split, 256):              %.4f ms (100 iters)\n", time_ms(start, end));
    
    start = clock();
    for (int i = 0; i < 100; i++) {
        BigUInt *c = biguint_mul_ntt_mont(a, b);
        biguint_free(c);
    }
    end = clock();
    printf("mul (NTT Montgomery, 256):         %.4f ms (100 iters)\n", time_ms(start, end));
    
    biguint_free(a);
    biguint_free(b);
    
    printf("\n================================\n");
    return 0;
}

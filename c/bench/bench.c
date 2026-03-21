#include <bigint.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#include "fft_split_bench.c"

int main() {
    printf("========== BENCHMARKS ==========\n\n");
    printf("ML-KEM / ML-DSA Vector Sizes\n");
    printf("=============================\n\n");
    
    srand(42);
    
    printf("ML-KEM-512 (256 words):\n");
    benchmark_fft_split("FFT", 256, 100);
    printf("\n");
    
    printf("ML-KEM-768 (512 words):\n");
    benchmark_fft_split("FFT", 512, 50);
    printf("\n");
    
    printf("ML-KEM-1024 (1024 words):\n");
    benchmark_fft_split("FFT", 1024, 20);
    printf("\n");
    
    printf("ML-DSA (2048 words):\n");
    benchmark_fft_split("FFT", 2048, 10);
    printf("\n");
    
    printf("Extended (3072 words):\n");
    benchmark_fft_split("FFT", 3072, 5);
    printf("\n");
    
    printf("Extended (4096 words):\n");
    benchmark_fft_split("FFT", 4096, 3);
    printf("\n");
    
    printf("================================\n");
    return 0;
}

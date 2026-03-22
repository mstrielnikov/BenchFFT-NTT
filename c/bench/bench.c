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

    benchmark_fft_split();

    printf("================================\n");
    return 0;
}

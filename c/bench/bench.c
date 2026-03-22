#include <bigint.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "fft_split_bench.h"
#include "fft_mersenne_bench.h"


int main() {
    printf("========== BENCHMARKS ==========\n\n");
    printf("ML-KEM / ML-DSA Vector Sizes\n");
    printf("=============================\n\n");

    srand(42);

    benchmark_fft_split();
    benchmark_fft_mersenne();

    printf("================================\n");
    return 0;
}

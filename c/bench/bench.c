#include <stdio.h>
#include "fft_split_bench.h"
#include "fft_mersenne_bench.h"
#include "ntt_mersenne_bench.h"
#include "ntt_mont.h"
#include "mult_standard_bench.h"


int main() {
    printf("========== BENCHMARKS ==========\n\n");
    printf("ML-KEM / ML-DSA Vector Sizes\n");
    printf("=============================\n\n");

    srand(42);

    benchmark_mult_standard();
    benchmark_fft_split();
    benchmark_fft_mersenne();
    benchmark_ntt_mersenne();
    benchmark_ntt_mont();

    printf("================================\n");
    return 0;
}

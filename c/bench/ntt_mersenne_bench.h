#pragma once
#include <bigint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <stdio.h>

#include "bench_utils.h"
#include "../src/ntt_mersenne.c"

#if HAS_AVX
#include "../src/ntt_mersenne_avx.c"
#endif


void benchmark_ntt_mersenne(void) {
    run_benchmark("NTT-M61", biguint_mul_ntt_mersenne);
#if HAS_AVX
    run_benchmark("NTT-M61-AVX", biguint_mul_ntt_mersenne_avx);
#endif
}

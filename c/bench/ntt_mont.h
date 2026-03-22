#pragma once
#include <bigint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <stdio.h>

#include "bench_utils.h"
#include "../src/ntt_mont.c"

#if HAS_AVX
#include "../src/ntt_mont_avx.c"
#endif


void benchmark_ntt_mont(void) {
    run_benchmark("NTT-Mont", biguint_mul_ntt_mont);
#if HAS_AVX
    run_benchmark("NTT-Mont-AVX", biguint_mul_ntt_mont_avx);
#endif
}

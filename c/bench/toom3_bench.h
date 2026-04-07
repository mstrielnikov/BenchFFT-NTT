#ifndef TOOM3_BENCH_H_INCLUDED
#define TOOM3_BENCH_H_INCLUDED

#include "bench_utils.h"
#include "../src/toom3.c"

#if HAS_AVX
#include "../src/toom3_avx.c"
#endif

void benchmark_toom3(void) {
    run_benchmark("Toom-Cook-3", biguint_mul_toom3);
#if HAS_AVX
    run_benchmark("Toom-Cook-3 AVX", biguint_mul_toom3_avx);
#endif
}

#endif /* TOOM3_BENCH_H_INCLUDED */

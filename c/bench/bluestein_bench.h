#ifndef BLUESTEIN_BENCH_H_INCLUDED
#define BLUESTEIN_BENCH_H_INCLUDED

#include "bench_utils.h"
#include "../src/bluestein.c"

void benchmark_bluestein(void) {
    run_benchmark("Bluestein CZT", biguint_mul_bluestein);
}

#endif /* BLUESTEIN_BENCH_H_INCLUDED */

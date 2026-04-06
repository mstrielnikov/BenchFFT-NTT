#ifndef NUSSBAUMER_BENCH_H_INCLUDED
#define NUSSBAUMER_BENCH_H_INCLUDED

#include "bench_utils.h"
#include "../src/nussbaumer.c"

void benchmark_nussbaumer(void) {
    run_benchmark("Nussbaumer", biguint_mul_nussbaumer);
}

#endif // NUSSBAUMER_BENCH_H_INCLUDED

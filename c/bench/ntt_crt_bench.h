#pragma once

#include "bench_utils.h"
#include "../src/ntt_crt.c"

void benchmark_ntt_crt(void) {
    run_benchmark("NTT (CRT)", biguint_mul_ntt_crt);
}

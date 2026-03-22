#include <bigint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <stdio.h>
#include "bench_utils.h"

#include "../src/fft_mersenne.c"


/* ── Per-algorithm wrappers (one per bench section) ─────────────────────── */
void benchmark_fft_mersenne(void) {
    run_benchmark("FFT-Mersenne", biguint_mul_fft_mersenne);
}

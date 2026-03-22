#include <bigint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <stdio.h>
#include "bench_utils.c"

#if HAS_AVX
#include "../src/fft_split_avx.c"
#endif


/* ── Per-algorithm wrappers (one per bench section) ─────────────────────── */
void benchmark_fft_split(void) {
    run_benchmark("FFT", biguint_mul_fft_split);

#if HAS_AVX
    run_benchmark("FFT-AVX", biguint_mul_fft_split_avx);
#endif
}

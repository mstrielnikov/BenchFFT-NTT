#include "bench_utils.h"

#include "../src/fft_mersenne.c"


/* ── Per-algorithm wrappers (one per bench section) ─────────────────────── */
void benchmark_fft_mersenne(void) {
    run_benchmark("FFT (M61)", biguint_mul_fft_mersenne);
}

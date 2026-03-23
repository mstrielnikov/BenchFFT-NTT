#include "bench_utils.h"
#include "../src/mult_standard.c"

#if HAS_AVX
#include "../src/mult_standard_avx.c"
#endif


/* ── Per-algorithm wrappers (one per bench section) ─────────────────────── */
void benchmark_mult_standard(void) {
    run_benchmark("STANDARD", biguint_mul_standard);

#if HAS_AVX
    run_benchmark("STANDARD AVX", biguint_mul_standard_avx);
#endif
}

#pragma once

#include "bench_utils.h"
#include "../src/ntt_mont.c"

#if HAS_AVX
#include "../src/ntt_mont_avx.c"
#endif

#ifdef __x86_64__
#include "../src/ntt_mont_asm.c"
#endif

void benchmark_ntt_mont(void) {
    run_benchmark("NTT (MONT)", biguint_mul_ntt_mont);
#if HAS_AVX
    run_benchmark("NTT (MONT) AVX", biguint_mul_ntt_mont_avx);
#endif
#ifdef __x86_64__
    run_benchmark("NTT (MONT) ASM", biguint_mul_ntt_mont_asm);
#endif
}

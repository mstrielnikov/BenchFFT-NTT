#pragma once
#include <bigint.h>
#include <stdlib.h>
#include <string.h>

#include "test_utils.h"
#include "../src/ntt_mont.c"
#include "../src/ntt_mont_m61.c"

#if HAS_AVX
#include "../src/ntt_mont_avx.c"
#endif

#ifdef __x86_64__
#include "../src/ntt_mont_asm.c"
#endif


void test_ntt_mont(void) {
    run_mul_tests("NTT (MONT)", biguint_mul_ntt_mont);
    run_mul_tests("NTT (MONT) M61", biguint_mul_ntt_mont_m61);
#if HAS_AVX
    run_mul_tests("NTT (MONT) AVX", biguint_mul_ntt_mont_avx);
#endif
#ifdef __x86_64__
    run_mul_tests("NTT (MONT) ASM", biguint_mul_ntt_mont_asm);
#endif
}

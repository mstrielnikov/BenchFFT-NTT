#pragma once
#include <bigint.h>
#include <stdlib.h>
#include <string.h>

#include "test_utils.h"
#include "../src/ntt_mersenne.c"

#if HAS_AVX
#include "../src/ntt_mersenne_avx.c"
#endif


void test_ntt_mersenne(void) {
    run_mul_tests("NTT MERSENNE", biguint_mul_ntt_mersenne);
#if HAS_AVX
    run_mul_tests("NTT MERSENNE AVX", biguint_mul_ntt_mersenne_avx);
#endif
}

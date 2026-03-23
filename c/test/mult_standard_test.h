#pragma once
#include <bigint.h>
#include <stdlib.h>
#include <string.h>

#include "test_utils.h"

#include "../src/mult_standard.c"

#if HAS_AVX
#include "../src/mult_standard_avx.c"
#endif

void test_mult_standard(void) {
    run_mul_tests("STANDARD", biguint_mul_standard);
#if HAS_AVX
    run_mul_tests("STANDARD AVX", biguint_mul_standard_avx);
#endif
}

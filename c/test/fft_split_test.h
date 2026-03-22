#include <bigint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "test_utils.h"
#include "../src/fft_split.c"

#if HAS_AVX
#include "../src/fft_split_avx.c"
#endif


void test_fft_split(void) {
    run_mul_tests("FFT SPLIT", biguint_mul_fft_split);
#if HAS_AVX
    run_mul_tests("FFT SPLIT AVX", biguint_mul_fft_split_avx);
#endif
}

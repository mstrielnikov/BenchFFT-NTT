#include <bigint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "test_utils.h"
#include "../src/fft_mersenne.c"


void test_fft_mersenne(void) {
    run_mul_tests("FFT MERSENNE", biguint_mul_fft_mersenne);
}

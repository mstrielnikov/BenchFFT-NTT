#ifndef TOOM3_TEST_H_INCLUDED
#define TOOM3_TEST_H_INCLUDED

#include "test_utils.h"
#include "../src/toom3.c"

#if HAS_AVX
#include "../src/toom3_avx.c"
#endif

void test_toom3(void) {
    run_mul_tests("Toom-Cook-3", biguint_mul_toom3);
#if HAS_AVX
    run_mul_tests("Toom-Cook-3 AVX", biguint_mul_toom3_avx);
#endif
}

#endif /* TOOM3_TEST_H_INCLUDED */

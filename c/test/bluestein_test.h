#ifndef BLUESTEIN_TEST_H_INCLUDED
#define BLUESTEIN_TEST_H_INCLUDED

#include "test_utils.h"
#include "../src/bluestein.c"

void test_bluestein(void) {
    run_mul_tests("Bluestein CZT", biguint_mul_bluestein);
}

#endif /* BLUESTEIN_TEST_H_INCLUDED */

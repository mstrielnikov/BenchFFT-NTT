#ifndef NUSSBAUMER_TEST_H_INCLUDED
#define NUSSBAUMER_TEST_H_INCLUDED

#include "test_utils.h"
#include "../src/nussbaumer.c"

void test_nussbaumer() {
    run_mul_tests("NUSSBAUMER FPT", biguint_mul_nussbaumer);
}

#endif // NUSSBAUMER_TEST_H_INCLUDED

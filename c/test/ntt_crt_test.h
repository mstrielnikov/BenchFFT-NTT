#pragma once
#include <bigint.h>
#include <stdlib.h>
#include <string.h>

#include "test_utils.h"
#include "../src/ntt_crt.c"

void test_ntt_crt(void) {
    run_mul_tests("NTT (CRT)", biguint_mul_ntt_crt);
}

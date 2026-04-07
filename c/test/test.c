#include <bigint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fft_split_test.h"
#include "fft_mersenne_test.h"
#include "ntt_mersenne_test.h"
#include "ntt_mont_test.h"
#include "mult_standard_test.h"
#include "ntt_crt_test.h"
#include "nussbaumer_test.h"
#include "bluestein_test.h"
#include "toom3_test.h"


int main() {
    test_mult_standard();
    test_fft_split();
    test_fft_mersenne();
    test_ntt_mersenne();
    test_ntt_mont();
    test_ntt_crt();
    test_nussbaumer();
    test_bluestein();
    test_toom3();
    
    printf("\n====================\n");
    printf("TOTAL: %d passed, %d failed\n", tests_passed, tests_failed);
    printf("====================\n");
    return tests_failed > 0 ? 1 : 0;
}

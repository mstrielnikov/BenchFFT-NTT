#include <bigint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fft_split_test.h"
#include "fft_mersenne_test.h"
#include "ntt_mersenne_test.h"


int main() {
    test_fft_split();
    test_fft_mersenne();
    test_ntt_mersenne();
    
    printf("\n====================\n");
    printf("TOTAL: %d passed, %d failed\n", tests_passed, tests_failed);
    printf("====================\n");
    return tests_failed > 0 ? 1 : 0;
}

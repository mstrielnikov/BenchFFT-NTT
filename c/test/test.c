#include <bigint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fft_split_test.c"

int main() {
    test_fft_split();
    
    printf("\n====================\n");
    printf("TOTAL: %d passed, %d failed\n", tests_passed, tests_failed);
    printf("====================\n");
    return tests_failed > 0 ? 1 : 0;
}

#include <stdio.h>
#include <stdint.h>
#include <math.h>

int main() {
    double max_val = (double)((uint64_t)-1LL) * (double)((uint64_t)-1LL) * 4096.0;
    printf("Max convolution value for 4096 words: %e\n", max_val);
    printf("Bits required: %f\n", log2(max_val));
    return 0;
}

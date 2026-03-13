#include <stdio.h>
#include <stdint.h>

int main() {
    uint64_t w1[4] = { -1ULL, -1ULL, -1ULL, -1ULL };
    printf("w1[0] = %lu\n", w1[0]);
    double d = (double)w1[0];
    uint64_t back = (uint64_t)d;
    printf("back = %lu\n", back);
    return 0;
}

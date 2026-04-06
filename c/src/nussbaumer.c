#ifndef NUSSBAUMER_C_INCLUDED
#define NUSSBAUMER_C_INCLUDED

#include "bigint.h"
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>


static int reverse_bits(int x, int p) {
    int rev = 0;
    for (int i = 0; i < p; i++) {
        rev = (rev << 1) | ((x >> i) & 1);
    }
    return rev;
}

static void nussbaumer_fft(__int128_t *A, int M, int K, bool inverse) {
    int p = 0;
    while ((1 << p) < M) p++;

    for (int i = 1; i < M; i++) {
        int j = reverse_bits(i, p);
        if (i < j) {
            for (int k = 0; k < K; k++) {
                __int128_t t = A[i * K + k];
                A[i * K + k] = A[j * K + k];
                A[j * K + k] = t;
            }
        }
    }

    int root_power = 2 * K / M;
    if (inverse) root_power = 2 * K - root_power;

    __int128_t *V = malloc(K * sizeof(__int128_t));

    for (int len = 2; len <= M; len <<= 1) {
        int half = len / 2;
        int step_power = (root_power * (M / len)) % (2 * K);

        for (int i = 0; i < M; i += len) {
            for (int j = 0; j < half; j++) {
                int shift = (j * step_power) % (2 * K);
                __int128_t *u_block = &A[(i + j) * K];
                __int128_t *v_block = &A[(i + j + half) * K];

                for (int k = 0; k < K; k++) {
                    int src = k - shift;
                    int negate = 0;
                    while (src < 0) { src += K; negate ^= 1; }
                    while (src >= K) { src -= K; negate ^= 1; }
                    V[k] = negate ? -v_block[src] : v_block[src];
                }

                for (int k = 0; k < K; k++) {
                    __int128_t u = u_block[k];
                    __int128_t v = V[k];
                    u_block[k] = u + v;
                    v_block[k] = u - v;
                }
            }
        }
    }

    if (inverse) {
        for (int i = 0; i < M * K; i++) {
            A[i] /= M;
        }
    }
    free(V);
}

static void poly_mul_mod(__int128_t *out, __int128_t *a, __int128_t *b, int K) {
    for (int i = 0; i < K; i++) out[i] = 0;
    for (int i = 0; i < K; i++) {
        if (a[i] == 0) continue;
        for (int j = 0; j < K; j++) {
            int pos = i + j;
            if (pos < K) {
                out[pos] += a[i] * b[j];
            } else {
                out[pos - K] -= a[i] * b[j];
            }
        }
    }
}

BigUInt *biguint_mul_nussbaumer(const BigUInt *a, const BigUInt *b) {
    if (biguint_is_zero(a) || biguint_is_zero(b)) return biguint_new();

    int La = a->len * 2;
    int Lb = b->len * 2;
    int L = 1;
    int M, K;
    
    // Find optimal L, K, M to prevent negacyclic wrap-around
    while (1) {
        int Ma = (La + L - 1) / L;
        int Mb = (Lb + L - 1) / L;
        int M_req = Ma + Mb - 1;
        if (M_req <= 0) M_req = 1;
        M = 1;
        while (M < M_req) M <<= 1;
        K = 2 * L;
        if (M <= 2 * K) break;
        L <<= 1;
    }

    __int128_t *A = calloc(M * K, sizeof(__int128_t));
    __int128_t *B = calloc(M * K, sizeof(__int128_t));

    for (int i = 0; i < La; i++) {
        uint32_t val = (i % 2 == 0) ? (uint32_t)(a->words[i / 2] & 0xFFFFFFFFULL) 
                                    : (uint32_t)(a->words[i / 2] >> 32);
        A[(i / L) * K + (i % L)] = val;
    }
    for (int i = 0; i < Lb; i++) {
        uint32_t val = (i % 2 == 0) ? (uint32_t)(b->words[i / 2] & 0xFFFFFFFFULL) 
                                    : (uint32_t)(b->words[i / 2] >> 32);
        B[(i / L) * K + (i % L)] = val;
    }

    nussbaumer_fft(A, M, K, false);
    nussbaumer_fft(B, M, K, false);

    __int128_t *C = calloc(M * K, sizeof(__int128_t));
    for (int m = 0; m < M; m++) {
        poly_mul_mod(&C[m * K], &A[m * K], &B[m * K], K);
    }

    nussbaumer_fft(C, M, K, true);

    // Overlap-add the blocks mapped back to linear 32-bit chunks
    size_t final_len32 = M * L + K;
    __int128_t *final_C = calloc(final_len32, sizeof(__int128_t));
    for (int m = 0; m < M; m++) {
        for (int k = 0; k < K; k++) {
            final_C[m * L + k] += C[m * K + k];
        }
    }

    BigUInt *res = biguint_new();
    size_t res_max_words = a->len + b->len;
    res->words = calloc(res_max_words, sizeof(uint64_t));
    res->len = res_max_words;
    res->cap = res_max_words;

    unsigned __int128 carry = 0;
    for (size_t i = 0; i < final_len32; i++) {
        unsigned __int128 val = (unsigned __int128)final_C[i] + carry;
        uint32_t chunk = (uint32_t)(val & 0xFFFFFFFFULL);
        carry = val >> 32;

        if (i / 2 < res->len) {
            if (i % 2 == 0) {
                res->words[i / 2] = chunk;
            } else {
                res->words[i / 2] |= ((uint64_t)chunk << 32);
            }
        }
    }
    
    size_t write_idx = final_len32 / 2;
    while (carry) {
        if (write_idx >= res->cap) {
            res->cap *= 2;
            res->words = realloc(res->words, res->cap * sizeof(uint64_t));
            memset(res->words + res->len, 0, (res->cap - res->len) * sizeof(uint64_t));
        }
        res->len = write_idx + 1;
        res->words[write_idx++] = (uint64_t)(carry & 0xFFFFFFFFFFFFFFFFULL);
        carry >>= 64;
    }

    free(A); free(B); free(C); free(final_C);

    biguint_normalize(res);
    return res;
}
#endif

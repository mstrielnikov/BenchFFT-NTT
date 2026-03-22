#include <bigint.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <m61_math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void fft_mersenne_inplace(double *re, double *im, size_t n, bool inverse) {
    for (size_t i = 1, j = 0; i < n; i++) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            double temp = re[i]; re[i] = re[j]; re[j] = temp;
            temp = im[i]; im[i] = im[j]; im[j] = temp;
        }
    }
    
    for (size_t len = 2; len <= n; len <<= 1) {
        double angle = (inverse ? -2.0 : 2.0) * M_PI / len;
        double wlen_re = cos(angle);
        double wlen_im = sin(angle);
        
        for (size_t i = 0; i < n; i += len) {
            double w_re = 1.0;
            double w_im = 0.0;
            
            for (size_t k = 0; k < len / 2; k++) {
                size_t idx = i + k;
                size_t idx2 = i + k + len / 2;
                
                double u_re = re[idx];
                double u_im = im[idx];
                double v_re = re[idx2] * w_re - im[idx2] * w_im;
                double v_im = re[idx2] * w_im + im[idx2] * w_re;
                
                re[idx] = u_re + v_re;
                im[idx] = u_im + v_im;
                re[idx2] = u_re - v_re;
                im[idx2] = u_im - v_im;
                
                double new_w_re = w_re * wlen_re - w_im * wlen_im;
                double new_w_im = w_re * wlen_im + w_im * wlen_re;
                w_re = new_w_re;
                w_im = new_w_im;
            }
        }
    }
    
    if (inverse) {
        for (size_t i = 0; i < n; i++) {
            re[i] /= n;
            im[i] /= n;
        }
    }
}

BigUInt *biguint_mul_fft_mersenne(const BigUInt *a, const BigUInt *b) {
    if (biguint_is_zero(a) || biguint_is_zero(b)) {
        return biguint_new();
    }
    
    size_t n = next_power_of_two(a->len + b->len - 1);
    
    double *re_a = calloc(n, sizeof(double));
    double *im_a = calloc(n, sizeof(double));
    double *re_b = calloc(n, sizeof(double));
    double *im_b = calloc(n, sizeof(double));
    
    if (!re_a || !im_a || !re_b || !im_b) {
        free(re_a); free(im_a); free(re_b); free(im_b);
        return NULL;
    }
    
    for (size_t i = 0; i < a->len; i++) re_a[i] = (double)a->words[i];
    for (size_t i = 0; i < b->len; i++) re_b[i] = (double)b->words[i];
    
    fft_mersenne_inplace(re_a, im_a, n, false);
    fft_mersenne_inplace(re_b, im_b, n, false);
    
    for (size_t i = 0; i < n; i++) {
        double r = re_a[i] * re_b[i] - im_a[i] * im_b[i];
        double im = re_a[i] * im_b[i] + im_a[i] * re_b[i];
        re_a[i] = r;
        im_a[i] = im;
    }
    
    fft_mersenne_inplace(re_a, im_a, n, true);
    
    BigUInt *res = biguint_new();
    if (!res) {
        free(re_a); free(im_a); free(re_b); free(im_b);
        return NULL;
    }
    
    for (size_t i = 0; i < n; i++) {
        re_a[i] += 0.5;
    }
    
    __uint128_t carry = 0;
    size_t len = 0;
    for (size_t i = 0; i < n; i++) {
        __uint128_t val = (__uint128_t)(re_a[i]) + carry;
        uint64_t word = (uint64_t)val;
        carry = val >> 64;
        
        uint64_t *tmp = realloc(res->words, (len + 1) * sizeof(uint64_t));
        if (!tmp) abort();
        res->words = tmp;
        res->words[len++] = m61_reduce(word);
    }
    
    while (carry > 0) {
        uint64_t word = (uint64_t)carry;
        uint64_t *tmp = realloc(res->words, (len + 1) * sizeof(uint64_t));
        if (!tmp) abort();
        res->words = tmp;
        res->words[len++] = m61_reduce(word);
        carry >>= M61_P;
    }
    
    res->len = len;
    res->cap = len;
    
    free(re_a); free(im_a); free(re_b); free(im_b);
    biguint_normalize(res);
    
    return res;
}

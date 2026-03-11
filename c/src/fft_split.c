#include <bigint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void fft_inplace(double *re, double *im, size_t n, bool inverse) {
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

static void *xrealloc(void *ptr, size_t size) {
    void *p = realloc(ptr, size);
    if (!p && size > 0) abort();
    return p;
}

BigUInt *biguint_mul_fft_split(const BigUInt *a, const BigUInt *b) {
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
    
    fft_inplace(re_a, im_a, n, false);
    fft_inplace(re_b, im_b, n, false);
    
    for (size_t i = 0; i < n; i++) {
        double r = re_a[i] * re_b[i] - im_a[i] * im_b[i];
        double im = re_a[i] * im_b[i] + im_a[i] * re_b[i];
        re_a[i] = r;
        im_a[i] = im;
    }
    
    fft_inplace(re_a, im_a, n, true);
    
    BigUInt *res = biguint_new();
    if (!res) {
        free(re_a); free(im_a); free(re_b); free(im_b);
        return NULL;
    }
    
    __uint128_t carry = 0;
    for (size_t i = 0; i < n; i++) {
        __uint128_t val = (__uint128_t)(re_a[i] + 0.5) + carry;
        res->words = xrealloc(res->words, (i + 1) * sizeof(uint64_t));
        res->words[i] = (uint64_t)val;
        carry = val >> 64;
    }
    
    while (carry > 0) {
        res->words = xrealloc(res->words, (res->len + 1) * sizeof(uint64_t));
        res->words[res->len++] = (uint64_t)carry;
        carry >>= 64;
    }
    
    res->len = n;
    res->cap = n;
    biguint_normalize(res);
    
    free(re_a); free(im_a); free(re_b); free(im_b);
    return res;
}

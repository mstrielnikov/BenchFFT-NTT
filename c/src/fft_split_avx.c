#include <bigint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <immintrin.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


#if HAS_AVX
static void fft_inplace_avx(double *re, double *im, size_t n, bool inverse) {
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
            
            size_t k = 0;
            if (len >= 8) {
                for (; k + 4 <= len / 2; k += 4) {
                    size_t idx = i + k;
                    
                    __m256d u_re = _mm256_loadu_pd(&re[idx]);
                    __m256d u_im = _mm256_loadu_pd(&im[idx]);
                    __m256d v_re = _mm256_loadu_pd(&re[idx + len / 2]);
                    __m256d v_im = _mm256_loadu_pd(&im[idx + len / 2]);
                    
                    __m256d v_re_new = _mm256_sub_pd(
                        _mm256_mul_pd(v_re, _mm256_set1_pd(w_re)),
                        _mm256_mul_pd(v_im, _mm256_set1_pd(w_im))
                    );
                    __m256d v_im_new = _mm256_add_pd(
                        _mm256_mul_pd(v_re, _mm256_set1_pd(w_im)),
                        _mm256_mul_pd(v_im, _mm256_set1_pd(w_re))
                    );
                    
                    __m256d sum_re = _mm256_add_pd(u_re, v_re_new);
                    __m256d sum_im = _mm256_add_pd(u_im, v_im_new);
                    __m256d diff_re = _mm256_sub_pd(u_re, v_re_new);
                    __m256d diff_im = _mm256_sub_pd(u_im, v_im_new);
                    
                    _mm256_storeu_pd(&re[idx], sum_re);
                    _mm256_storeu_pd(&im[idx], sum_im);
                    _mm256_storeu_pd(&re[idx + len / 2], diff_re);
                    _mm256_storeu_pd(&im[idx + len / 2], diff_im);
                    
                    double new_w_re = w_re * wlen_re - w_im * wlen_im;
                    double new_w_im = w_re * wlen_im + w_im * wlen_re;
                    w_re = new_w_re;
                    w_im = new_w_im;
                }
            }
            
            for (; k < len / 2; k++) {
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
        double inv_n = 1.0 / n;
        __m256d inv_n_vec = _mm256_set1_pd(inv_n);
        for (size_t i = 0; i < n; i += 4) {
            size_t remain = n - i;
            if (remain >= 4) {
                __m256d re_vec = _mm256_loadu_pd(&re[i]);
                __m256d im_vec = _mm256_loadu_pd(&im[i]);
                re_vec = _mm256_mul_pd(re_vec, inv_n_vec);
                im_vec = _mm256_mul_pd(im_vec, inv_n_vec);
                _mm256_storeu_pd(&re[i], re_vec);
                _mm256_storeu_pd(&im[i], im_vec);
            } else {
                for (size_t j = i; j < n; j++) {
                    re[j] *= inv_n;
                    im[j] *= inv_n;
                }
                break;
            }
        }
    }
}

BigUInt *biguint_mul_fft_split_avx(const BigUInt *a, const BigUInt *b) {
    if (biguint_is_zero(a) || biguint_is_zero(b)) {
        return biguint_new();
    }
    
    size_t n = next_power_of_two(a->len + b->len - 1);
    size_t result_len = a->len + b->len - 1;
    
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
    
    fft_inplace_avx(re_a, im_a, n, false);
    fft_inplace_avx(re_b, im_b, n, false);
    
    for (size_t i = 0; i < n; i += 4) {
        size_t remain = n - i;
        if (remain >= 4) {
            __m256d re_a_vec = _mm256_loadu_pd(&re_a[i]);
            __m256d im_a_vec = _mm256_loadu_pd(&im_a[i]);
            __m256d re_b_vec = _mm256_loadu_pd(&re_b[i]);
            __m256d im_b_vec = _mm256_loadu_pd(&im_b[i]);
            
            __m256d r_re = _mm256_sub_pd(
                _mm256_mul_pd(re_a_vec, re_b_vec),
                _mm256_mul_pd(im_a_vec, im_b_vec)
            );
            __m256d r_im = _mm256_add_pd(
                _mm256_mul_pd(re_a_vec, im_b_vec),
                _mm256_mul_pd(im_a_vec, re_b_vec)
            );
            
            _mm256_storeu_pd(&re_a[i], r_re);
            _mm256_storeu_pd(&im_a[i], r_im);
        } else {
            for (size_t j = i; j < n; j++) {
                double r = re_a[j] * re_b[j] - im_a[j] * im_b[j];
                double im = re_a[j] * im_b[j] + im_a[j] * re_b[j];
                re_a[j] = r;
                im_a[j] = im;
            }
        }
    }
    
    fft_inplace_avx(re_a, im_a, n, true);
    
    BigUInt *res = biguint_new();
    if (!res) {
        free(re_a); free(im_a); free(re_b); free(im_b);
        return NULL;
    }
    
    res->words = calloc(result_len, sizeof(uint64_t));
    if (!res->words) {
        free(re_a); free(im_a); free(re_b); free(im_b);
        free(res);
        return NULL;
    }
    res->len = result_len;
    res->cap = result_len;
    
    __uint128_t carry = 0;
    for (size_t i = 0; i < n; i++) {
        __uint128_t val = (__uint128_t)(re_a[i] + 0.5) + carry;
        if (i < result_len) {
            res->words[i] = (uint64_t)val;
        }
        carry = val >> 64;
    }
    
    free(re_a); free(im_a); free(re_b); free(im_b);
    
    biguint_normalize(res);
    return res;
}
#endif

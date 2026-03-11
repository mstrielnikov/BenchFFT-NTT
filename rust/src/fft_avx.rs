use crate::bigint::{next_power_of_two, BigUInt};
use std::f64::consts::PI;

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

#[cfg(target_arch = "x86_64")]
unsafe fn fft_inplace_avx(re: &mut [f64], im: &mut [f64], inverse: bool) {
    let n = re.len();

    let mut j = 0;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            re.swap(i, j);
            im.swap(i, j);
        }
    }

    let mut len = 2;
    while len <= n {
        let angle = if inverse { -2.0 } else { 2.0 } * PI / len as f64;
        let wlen_re = angle.cos();
        let wlen_im = angle.sin();

        for i in (0..n).step_by(len) {
            let mut w_re = 1.0;
            let mut w_im = 0.0;

            let mut k = 0;
            if len >= 8 {
                for k_outer in (0..len / 2).step_by(4) {
                    let idx = i + k_outer;
                    let idx2 = i + k_outer + len / 2;

                    let u_re = _mm256_loadu_pd(&re[idx]);
                    let u_im = _mm256_loadu_pd(&im[idx]);
                    let v_re = _mm256_loadu_pd(&re[idx2]);
                    let v_im = _mm256_loadu_pd(&im[idx2]);

                    let w_re_vec = _mm256_set1_pd(w_re);
                    let w_im_vec = _mm256_set1_pd(w_im);

                    let vr_wr = _mm256_mul_pd(v_re, w_re_vec);
                    let vi_wi = _mm256_mul_pd(v_im, w_im_vec);
                    let vr_wi = _mm256_mul_pd(v_re, w_im_vec);
                    let vi_wr = _mm256_mul_pd(v_im, w_re_vec);

                    let v_new_re = _mm256_sub_pd(vr_wr, vi_wi);
                    let v_new_im = _mm256_add_pd(vr_wi, vi_wr);

                    let sum_re = _mm256_add_pd(u_re, v_new_re);
                    let sum_im = _mm256_add_pd(u_im, v_new_im);
                    let diff_re = _mm256_sub_pd(u_re, v_new_re);
                    let diff_im = _mm256_sub_pd(u_im, v_new_im);

                    _mm256_storeu_pd(&mut re[idx], sum_re);
                    _mm256_storeu_pd(&mut im[idx], sum_im);
                    _mm256_storeu_pd(&mut re[idx2], diff_re);
                    _mm256_storeu_pd(&mut im[idx2], diff_im);

                    let new_w_re = w_re * wlen_re - w_im * wlen_im;
                    let new_w_im = w_re * wlen_im + w_im * wlen_re;
                    w_re = new_w_re;
                    w_im = new_w_im;
                }
                k = (len / 2 / 4) * 4;
            }

            for k_outer in k..len / 2 {
                let idx = i + k_outer;
                let idx2 = i + k_outer + len / 2;

                let u_re = re[idx];
                let u_im = im[idx];
                let v_re = re[idx2] * w_re - im[idx2] * w_im;
                let v_im = re[idx2] * w_im + im[idx2] * w_re;

                re[idx] = u_re + v_re;
                im[idx] = u_im + v_im;
                re[idx2] = u_re - v_re;
                im[idx2] = u_im - v_im;

                let new_w_re = w_re * wlen_re - w_im * wlen_im;
                let new_w_im = w_re * wlen_im + w_im * wlen_re;
                w_re = new_w_re;
                w_im = new_w_im;
            }
        }
        len <<= 1;
    }

    if inverse {
        let inv_n = 1.0 / n as f64;
        let inv_n_vec = _mm256_set1_pd(inv_n);
        for i in (0..n).step_by(4) {
            if i + 4 <= n {
                let re_vec = _mm256_loadu_pd(&re[i]);
                let im_vec = _mm256_loadu_pd(&im[i]);
                let re_scaled = _mm256_mul_pd(re_vec, inv_n_vec);
                let im_scaled = _mm256_mul_pd(im_vec, inv_n_vec);
                _mm256_storeu_pd(&mut re[i], re_scaled);
                _mm256_storeu_pd(&mut im[i], im_scaled);
            } else {
                for j in i..n {
                    re[j] *= inv_n;
                    im[j] *= inv_n;
                }
            }
        }
    }
}

#[cfg(target_arch = "x86_64")]
pub fn biguint_mul_fft_split_avx(a: &BigUInt, b: &BigUInt) -> BigUInt {
    if a.is_zero() || b.is_zero() {
        return BigUInt::new();
    }

    let n = next_power_of_two(a.words().len() + b.words().len() - 1);

    let mut re_a = vec![0.0; n];
    let mut im_a = vec![0.0; n];
    let mut re_b = vec![0.0; n];
    let mut im_b = vec![0.0; n];

    for i in 0..a.words().len() {
        re_a[i] = a.words()[i] as f64;
    }
    for i in 0..b.words().len() {
        re_b[i] = b.words()[i] as f64;
    }

    unsafe {
        fft_inplace_avx(&mut re_a, &mut im_a, false);
        fft_inplace_avx(&mut re_b, &mut im_b, false);
    }

    unsafe {
        for i in (0..n).step_by(4) {
            if i + 4 <= n {
                let re_a_vec = _mm256_loadu_pd(&re_a[i]);
                let im_a_vec = _mm256_loadu_pd(&im_a[i]);
                let re_b_vec = _mm256_loadu_pd(&re_b[i]);
                let im_b_vec = _mm256_loadu_pd(&im_b[i]);

                let vr_wr = _mm256_mul_pd(re_a_vec, re_b_vec);
                let vi_wi = _mm256_mul_pd(im_a_vec, im_b_vec);
                let r_re = _mm256_sub_pd(vr_wr, vi_wi);

                let vr_wi = _mm256_mul_pd(re_a_vec, im_b_vec);
                let vi_wr = _mm256_mul_pd(im_a_vec, re_b_vec);
                let r_im = _mm256_add_pd(vr_wi, vi_wr);

                _mm256_storeu_pd(&mut re_a[i], r_re);
                _mm256_storeu_pd(&mut im_a[i], r_im);
            } else {
                for j in i..n {
                    let r = re_a[j] * re_b[j] - im_a[j] * im_b[j];
                    let im = re_a[j] * im_b[j] + im_a[j] * re_b[j];
                    re_a[j] = r;
                    im_a[j] = im;
                }
            }
        }
    }

    unsafe {
        fft_inplace_avx(&mut re_a, &mut im_a, true);
    }

    let mut result = Vec::with_capacity(n);
    let mut carry: u128 = 0;
    for i in 0..n {
        let val = (re_a[i] + 0.5) as u128 + carry;
        result.push(val as u64);
        carry = val >> 64;
    }

    while carry > 0 {
        result.push(carry as u64);
        carry >>= 64;
    }

    BigUInt::from_slice(&result)
}

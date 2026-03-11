use std::f64::consts::PI;

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

#[derive(Clone, Debug, PartialEq)]
pub struct BigUInt {
    words: Vec<u64>,
}

impl BigUInt {
    pub fn new() -> Self {
        BigUInt { words: Vec::new() }
    }

    pub fn from_u64(n: u64) -> Self {
        if n == 0 {
            BigUInt::new()
        } else {
            BigUInt { words: vec![n] }
        }
    }

    pub fn from_slice(words: &[u64]) -> Self {
        let mut result = BigUInt {
            words: words.to_vec(),
        };
        result.normalize();
        result
    }

    pub fn is_zero(&self) -> bool {
        self.words.is_empty()
    }

    fn normalize(&mut self) {
        while let Some(&last) = self.words.last() {
            if last == 0 {
                self.words.pop();
            } else {
                break;
            }
        }
    }

    pub fn len(&self) -> usize {
        self.words.len()
    }

    pub fn words(&self) -> &[u64] {
        &self.words
    }
}

fn next_power_of_two(n: usize) -> usize {
    let mut n = n;
    n -= 1;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n |= n >> 32;
    n + 1
}

fn fft_inplace(re: &mut [f64], im: &mut [f64], inverse: bool) {
    let n = re.len();

    // Bit-reversal permutation
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

    // Cooley-Tukey FFT
    let mut len = 2;
    while len <= n {
        let angle = if inverse { -2.0 } else { 2.0 } * PI / len as f64;
        let wlen_re = angle.cos();
        let wlen_im = angle.sin();

        for i in (0..n).step_by(len) {
            let mut w_re = 1.0;
            let mut w_im = 0.0;

            for k in 0..len / 2 {
                let idx = i + k;
                let idx2 = i + k + len / 2;

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
        let n = n as f64;
        for i in 0..n as usize {
            re[i] /= n;
            im[i] /= n;
        }
    }
}

#[cfg(target_arch = "x86_64")]
unsafe fn fft_inplace_avx(re: &mut [f64], im: &mut [f64], inverse: bool) {
    let n = re.len();

    // Bit-reversal permutation
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

    // Cooley-Tukey FFT with AVX vectorization
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

                    // Load 4 doubles (256 bits = 4 f64)
                    let u_re = _mm256_loadu_pd(&re[idx]);
                    let u_im = _mm256_loadu_pd(&im[idx]);
                    let v_re = _mm256_loadu_pd(&re[idx2]);
                    let v_im = _mm256_loadu_pd(&im[idx2]);

                    // Multiply v by w (complex multiplication)
                    // v * w = (v_re * w_re - v_im * w_im) + i(v_re * w_im + v_im * w_re)
                    let w_re_vec = _mm256_set1_pd(w_re);
                    let w_im_vec = _mm256_set1_pd(w_im);

                    // v_re * w_re
                    let vr_wr = _mm256_mul_pd(v_re, w_re_vec);
                    // v_im * w_im
                    let vi_wi = _mm256_mul_pd(v_im, w_im_vec);
                    // v_re * w_im
                    let vr_wi = _mm256_mul_pd(v_re, w_im_vec);
                    // v_im * w_re
                    let vi_wr = _mm256_mul_pd(v_im, w_re_vec);

                    // v_new_re = v_re * w_re - v_im * w_im
                    let v_new_re = _mm256_sub_pd(vr_wr, vi_wi);
                    // v_new_im = v_re * w_im + v_im * w_re
                    let v_new_im = _mm256_add_pd(vr_wi, vi_wr);

                    // u + v
                    let sum_re = _mm256_add_pd(u_re, v_new_re);
                    let sum_im = _mm256_add_pd(u_im, v_new_im);
                    // u - v
                    let diff_re = _mm256_sub_pd(u_re, v_new_re);
                    let diff_im = _mm256_sub_pd(u_im, v_new_im);

                    _mm256_storeu_pd(&mut re[idx], sum_re);
                    _mm256_storeu_pd(&mut im[idx], sum_im);
                    _mm256_storeu_pd(&mut re[idx2], diff_re);
                    _mm256_storeu_pd(&mut im[idx2], diff_im);

                    // Update w for next 4 elements
                    let new_w_re = w_re * wlen_re - w_im * wlen_im;
                    let new_w_im = w_re * wlen_im + w_im * wlen_re;
                    w_re = new_w_re;
                    w_im = new_w_im;
                }
                k = (len / 2 / 4) * 4;
            }

            // Scalar tail
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

            // Scalar tail
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

pub fn biguint_mul_fft_split(a: &BigUInt, b: &BigUInt) -> BigUInt {
    if a.is_zero() || b.is_zero() {
        return BigUInt::new();
    }

    let n = next_power_of_two(a.words.len() + b.words.len() - 1);

    let mut re_a = vec![0.0; n];
    let mut im_a = vec![0.0; n];
    let mut re_b = vec![0.0; n];
    let mut im_b = vec![0.0; n];

    for i in 0..a.words.len() {
        re_a[i] = a.words[i] as f64;
    }
    for i in 0..b.words.len() {
        re_b[i] = b.words[i] as f64;
    }

    fft_inplace(&mut re_a, &mut im_a, false);
    fft_inplace(&mut re_b, &mut im_b, false);

    for i in 0..n {
        let r = re_a[i] * re_b[i] - im_a[i] * im_b[i];
        let im = re_a[i] * im_b[i] + im_a[i] * re_b[i];
        re_a[i] = r;
        im_a[i] = im;
    }

    fft_inplace(&mut re_a, &mut im_a, true);

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

#[cfg(target_arch = "x86_64")]
pub fn biguint_mul_fft_split_avx(a: &BigUInt, b: &BigUInt) -> BigUInt {
    if a.is_zero() || b.is_zero() {
        return BigUInt::new();
    }

    let n = next_power_of_two(a.words.len() + b.words.len() - 1);

    let mut re_a = vec![0.0; n];
    let mut im_a = vec![0.0; n];
    let mut re_b = vec![0.0; n];
    let mut im_b = vec![0.0; n];

    for i in 0..a.words.len() {
        re_a[i] = a.words[i] as f64;
    }
    for i in 0..b.words.len() {
        re_b[i] = b.words[i] as f64;
    }

    unsafe {
        fft_inplace_avx(&mut re_a, &mut im_a, false);
        fft_inplace_avx(&mut re_b, &mut im_b, false);
    }

    // AVX pointwise multiply
    unsafe {
        for i in (0..n).step_by(4) {
            if i + 4 <= n {
                let re_a_vec = _mm256_loadu_pd(&re_a[i]);
                let im_a_vec = _mm256_loadu_pd(&im_a[i]);
                let re_b_vec = _mm256_loadu_pd(&re_b[i]);
                let im_b_vec = _mm256_loadu_pd(&im_b[i]);

                // r_re = re_a * re_b - im_a * im_b
                let vr_wr = _mm256_mul_pd(re_a_vec, re_b_vec);
                let vi_wi = _mm256_mul_pd(im_a_vec, im_b_vec);
                let r_re = _mm256_sub_pd(vr_wr, vi_wi);

                // r_im = re_a * im_b + im_a * re_b
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

// Montgomery NTT
const MONT_MOD: u64 = 998244353;
const MONT_ROOT: u64 = 3;

fn mod_pow(mut a: u64, mut e: u64, m: u64) -> u64 {
    let mut res = 1;
    while e > 0 {
        if e & 1 == 1 {
            res = ((res as u128 * a as u128) % m as u128) as u64;
        }
        a = ((a as u128 * a as u128) % m as u128) as u64;
        e >>= 1;
    }
    res
}

fn mod_mul(a: u64, b: u64) -> u64 {
    ((a as u128 * b as u128) % MONT_MOD as u128) as u64
}

fn ntt_inplace(a: &mut [u64], inverse: bool) {
    let n = a.len();

    // Bit-reversal
    let mut j = 0;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            a.swap(i, j);
        }
    }

    let mut len = 2;
    while len <= n {
        let wlen = mod_pow(MONT_ROOT, (MONT_MOD - 1) / len as u64, MONT_MOD);
        let wlen = if inverse {
            mod_pow(wlen, MONT_MOD - 2, MONT_MOD)
        } else {
            wlen
        };

        for i in (0..n).step_by(len) {
            let mut w = 1u64;
            for j in 0..len / 2 {
                let u = a[i + j];
                let v = mod_mul(a[i + j + len / 2], w);
                a[i + j] = u + v;
                if a[i + j] >= MONT_MOD {
                    a[i + j] -= MONT_MOD;
                }
                a[i + j + len / 2] = u + MONT_MOD - v;
                if a[i + j + len / 2] >= MONT_MOD {
                    a[i + j + len / 2] -= MONT_MOD;
                }
                w = mod_mul(w, wlen);
            }
        }
        len <<= 1;
    }

    if inverse {
        let inv_n = mod_pow(n as u64, MONT_MOD - 2, MONT_MOD);
        for x in a.iter_mut() {
            *x = mod_mul(*x, inv_n);
        }
    }
}

pub fn biguint_mul_ntt_mont(a: &BigUInt, b: &BigUInt) -> BigUInt {
    if a.is_zero() || b.is_zero() {
        return BigUInt::new();
    }

    // For NTT, we need to work modulo MONT_MOD
    // Convert words to mod MONT_MOD representation
    let n = next_power_of_two(a.words.len() + b.words.len() - 1);

    let mut fa: Vec<u64> = vec![0; n];
    let mut fb: Vec<u64> = vec![0; n];

    for i in 0..a.words.len() {
        fa[i] = a.words[i] % MONT_MOD;
    }
    for i in 0..b.words.len() {
        fb[i] = b.words[i] % MONT_MOD;
    }

    // NTT
    ntt_inplace(&mut fa, false);
    ntt_inplace(&mut fb, false);

    // Pointwise multiplication
    for i in 0..n {
        fa[i] = mod_mul(fa[i], fb[i]);
    }

    // Inverse NTT
    ntt_inplace(&mut fa, true);

    // Convert back to BigUInt (words are already in range 0..MONT_MOD)
    let result = BigUInt::from_slice(&fa);
    result
}

#[cfg(target_arch = "x86_64")]
pub fn biguint_mul_ntt_mont_avx(a: &BigUInt, b: &BigUInt) -> BigUInt {
    if a.is_zero() || b.is_zero() {
        return BigUInt::new();
    }

    let n = next_power_of_two(a.words.len() + b.words.len() - 1);

    let mut fa: Vec<u64> = vec![0; n];
    let mut fb: Vec<u64> = vec![0; n];

    for i in 0..a.words.len() {
        fa[i] = a.words[i] % MONT_MOD;
    }
    for i in 0..b.words.len() {
        fb[i] = b.words[i] % MONT_MOD;
    }

    ntt_inplace(&mut fa, false);
    ntt_inplace(&mut fb, false);

    for i in 0..n {
        fa[i] = mod_mul(fa[i], fb[i]);
    }

    ntt_inplace(&mut fa, true);

    let result = BigUInt::from_slice(&fa);
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mul_small() {
        let a = BigUInt::from_u64(12345);
        let b = BigUInt::from_u64(67890);
        let c = biguint_mul_fft_split(&a, &b);
        let expected = BigUInt::from_u64(12345 * 67890);
        assert_eq!(c, expected);
    }

    #[test]
    fn test_mul_large() {
        let words1: Vec<u64> = (0..16).map(|i| i as u64).collect();
        let words2: Vec<u64> = (0..16).map(|i| (i * 2) as u64).collect();
        let a = BigUInt::from_slice(&words1);
        let b = BigUInt::from_slice(&words2);
        let _c = biguint_mul_fft_split(&a, &b);
    }

    #[test]
    fn test_ntt_small() {
        let a = BigUInt::from_u64(12345);
        let b = BigUInt::from_u64(67890);
        let c = biguint_mul_ntt_mont(&a, &b);
        let expected = BigUInt::from_u64(12345 * 67890);
        assert_eq!(c, expected);
    }
}

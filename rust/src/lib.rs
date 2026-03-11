use std::f64::consts::PI;

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

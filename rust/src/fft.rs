use crate::bigint::{next_power_of_two, BigUInt};
use std::f64::consts::PI;

fn fft_inplace(re: &mut [f64], im: &mut [f64], inverse: bool) {
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
}

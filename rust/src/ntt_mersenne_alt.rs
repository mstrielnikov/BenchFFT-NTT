use crate::bigint::{next_power_of_two, BigUInt};

const M61_MOD: u64 = 2305843009213693951;
const M61_P: u64 = 61;
const NTT_MOD: u64 = 998244353;
const NTT_ROOT: u64 = 3;

#[inline]
fn m61_reduce(x: u64) -> u64 {
    let t = (x >> M61_P) + (x & M61_MOD);
    let t = t + (t >> M61_P);
    t & M61_MOD
}

fn mod_pow(mut a: u64, mut e: u64, m: u64) -> u64 {
    let mut res = 1;
    while e > 0 {
        if e & 1 == 1 {
            res = ((res as u128) * (a as u128) % m as u128) as u64;
        }
        a = ((a as u128) * (a as u128) % m as u128) as u64;
        e >>= 1;
    }
    res
}

#[inline]
fn mod_mul(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128) * (b as u128) % m as u128) as u64
}

fn ntt_inplace(a: &mut [u64], mod_val: u64, root: u64, inverse: bool) {
    let n = a.len();

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
        let wlen = mod_pow(root, (mod_val - 1) / len as u64, mod_val);
        let wlen = if inverse {
            mod_pow(wlen, mod_val - 2, mod_val)
        } else {
            wlen
        };

        for i in (0..n).step_by(len) {
            let mut w = 1u64;
            for j in 0..len / 2 {
                let u = a[i + j];
                let v = mod_mul(a[i + j + len / 2], w, mod_val);
                a[i + j] = u + v;
                if a[i + j] >= mod_val {
                    a[i + j] -= mod_val;
                }
                a[i + j + len / 2] = u + mod_val - v;
                if a[i + j + len / 2] >= mod_val {
                    a[i + j + len / 2] -= mod_val;
                }
                w = mod_mul(w, wlen, mod_val);
            }
        }
        len <<= 1;
    }

    if inverse {
        let inv_n = mod_pow(n as u64, mod_val - 2, mod_val);
        for x in a.iter_mut() {
            *x = mod_mul(*x, inv_n, mod_val);
        }
    }
}

pub fn biguint_mul_ntt_mersenne_alt(a: &BigUInt, b: &BigUInt) -> BigUInt {
    if a.is_zero() || b.is_zero() {
        return BigUInt::new();
    }

    let n = next_power_of_two(a.words().len() + b.words().len() - 1);
    let result_len = a.words().len() + b.words().len() - 1;

    let mut fa: Vec<u64> = vec![0; n];
    let mut fb: Vec<u64> = vec![0; n];

    for i in 0..a.words().len() {
        fa[i] = a.words()[i] % NTT_MOD;
    }
    for i in 0..b.words().len() {
        fb[i] = b.words()[i] % NTT_MOD;
    }

    ntt_inplace(&mut fa, NTT_MOD, NTT_ROOT, false);
    ntt_inplace(&mut fb, NTT_MOD, NTT_ROOT, false);

    for i in 0..n {
        fa[i] = mod_mul(fa[i], fb[i], NTT_MOD);
    }

    ntt_inplace(&mut fa, NTT_MOD, NTT_ROOT, true);

    let mut result = Vec::with_capacity(result_len);
    for i in 0..result_len {
        result.push(m61_reduce(fa[i]));
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
        let c = biguint_mul_ntt_mersenne_alt(&a, &b);
        let expected = BigUInt::from_u64(12345 * 67890);
        assert_eq!(c, expected);
    }
}

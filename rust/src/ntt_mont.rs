use crate::bigint::{next_power_of_two, BigUInt};

const MONT_MOD: u64 = 998244353;
const MONT_ROOT: u64 = 3;

const MONTGOMERY_R_INV: u64 = 17450252288407896063; // (-998244353)^-1 mod 2^64
const R2_MOD: u64 = 299560064;

#[inline]
fn montgomery_mul(a: u64, b: u64) -> u64 {
    let t = (a as u128) * (b as u128);
    let m = (t as u64).wrapping_mul(MONTGOMERY_R_INV);
    let m_mod = (m as u128) * (MONT_MOD as u128);
    let mut result = ((t + m_mod) >> 64) as u64;
    if result >= MONT_MOD {
        result -= MONT_MOD;
    }
    result
}

fn mod_pow_mont(mut base: u64, mut exp: u64) -> u64 {
    let mut result = montgomery_mul(1, R2_MOD);
    while exp > 0 {
        if exp & 1 == 1 {
            result = montgomery_mul(result, base);
        }
        base = montgomery_mul(base, base);
        exp >>= 1;
    }
    result
}

fn ntt_inplace(a: &mut [u64], root: u64, inverse: bool) {
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
        let wlen = mod_pow_mont(root, (MONT_MOD - 1) / len as u64);
        let wlen = if inverse {
            mod_pow_mont(wlen, MONT_MOD - 2)
        } else {
            wlen
        };

        for i in (0..n).step_by(len) {
            let mut w = montgomery_mul(1, R2_MOD);
            for j in 0..len / 2 {
                let u = a[i + j];
                let v = montgomery_mul(a[i + j + len / 2], w);
                a[i + j] = u + v;
                if a[i + j] >= MONT_MOD {
                    a[i + j] -= MONT_MOD;
                }
                a[i + j + len / 2] = u + MONT_MOD - v;
                if a[i + j + len / 2] >= MONT_MOD {
                    a[i + j + len / 2] -= MONT_MOD;
                }
                w = montgomery_mul(w, wlen);
            }
        }
        len <<= 1;
    }

    if inverse {
        let inv_n = mod_pow_mont(montgomery_mul(n as u64, R2_MOD), MONT_MOD - 2);
        for x in a.iter_mut() {
            *x = montgomery_mul(*x, inv_n);
        }
    }
}

pub fn biguint_mul_ntt_mont(a: &BigUInt, b: &BigUInt) -> BigUInt {
    if a.is_zero() || b.is_zero() {
        return BigUInt::new();
    }

    let n = next_power_of_two(a.words().len() + b.words().len() - 1);

    let mut fa: Vec<u64> = vec![0; n];
    let mut fb: Vec<u64> = vec![0; n];

    for i in 0..a.words().len() {
        fa[i] = montgomery_mul(a.words()[i] % MONT_MOD, R2_MOD);
    }
    for i in 0..b.words().len() {
        fb[i] = montgomery_mul(b.words()[i] % MONT_MOD, R2_MOD);
    }

    let mont_root = montgomery_mul(MONT_ROOT, R2_MOD);

    ntt_inplace(&mut fa, mont_root, false);
    ntt_inplace(&mut fb, mont_root, false);

    for i in 0..n {
        fa[i] = montgomery_mul(fa[i], fb[i]);
    }

    ntt_inplace(&mut fa, mont_root, true);

    let mut result = Vec::with_capacity(n);
    for i in 0..n {
        result.push(montgomery_mul(fa[i], 1));
    }

    BigUInt::from_slice(&result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ntt_small() {
        let a = BigUInt::from_u64(12345);
        let b = BigUInt::from_u64(67890);
        let c = biguint_mul_ntt_mont(&a, &b);
        let expected = BigUInt::from_u64(12345 * 67890);
        assert_eq!(c, expected);
    }
}

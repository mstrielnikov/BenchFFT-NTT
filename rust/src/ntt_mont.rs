use crate::bigint::{next_power_of_two, BigUInt};

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

#[inline]
fn mod_mul(a: u64, b: u64) -> u64 {
    ((a as u128 * b as u128) % MONT_MOD as u128) as u64
}

fn ntt_inplace(a: &mut [u64], inverse: bool) {
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

    let n = next_power_of_two(a.words().len() + b.words().len() - 1);

    let mut fa: Vec<u64> = vec![0; n];
    let mut fb: Vec<u64> = vec![0; n];

    for i in 0..a.words().len() {
        fa[i] = a.words()[i] % MONT_MOD;
    }
    for i in 0..b.words().len() {
        fb[i] = b.words()[i] % MONT_MOD;
    }

    ntt_inplace(&mut fa, false);
    ntt_inplace(&mut fb, false);

    for i in 0..n {
        fa[i] = mod_mul(fa[i], fb[i]);
    }

    ntt_inplace(&mut fa, true);

    BigUInt::from_slice(&fa)
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

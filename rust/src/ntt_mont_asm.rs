use crate::bigint::{next_power_of_two, BigUInt};

const NTT_MOD: u64 = 998244353;
const NTT_ROOT: u64 = 3;
const MONTGOMERY_R_INV: u64 = 330382944;

#[inline]
fn montgomery_mul(a: u64, b: u64) -> u64 {
    let lo = (a as u128) * (b as u128);
    let mi = ((lo & 0xffffffff) as u64).wrapping_mul(MONTGOMERY_R_INV);
    let hi = ((lo >> 64) as u64).wrapping_add(((mi as u128 * NTT_MOD as u128) >> 64) as u64);
    let result = hi.wrapping_sub(NTT_MOD);
    if result > hi {
        result.wrapping_add(NTT_MOD)
    } else {
        result
    }
}

#[inline]
fn mod_mul(a: u64, b: u64) -> u64 {
    ((a as u128) * (b as u128) % NTT_MOD as u128) as u64
}

fn mod_pow_asm(mut base: u64, mut exp: u64) -> u64 {
    let mut result = 1u64;
    base %= NTT_MOD;

    while exp > 0 {
        if exp & 1 == 1 {
            result = mod_mul(result, base);
        }
        base = mod_mul(base, base);
        exp >>= 1;
    }
    result
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
        let wlen = mod_pow_asm(root, (mod_val - 1) / len as u64);
        let wlen = if inverse {
            mod_pow_asm(wlen, mod_val - 2)
        } else {
            wlen
        };

        for i in (0..n).step_by(len) {
            let mut w = 1u64;
            for j in 0..len / 2 {
                let u = a[i + j];
                let v = mod_mul(a[i + j + len / 2], w);
                a[i + j] = u + v;
                if a[i + j] >= mod_val {
                    a[i + j] -= mod_val;
                }
                a[i + j + len / 2] = u + mod_val - v;
                if a[i + j + len / 2] >= mod_val {
                    a[i + j + len / 2] -= mod_val;
                }
                w = mod_mul(w, wlen);
            }
        }
        len <<= 1;
    }

    if inverse {
        let inv_n = mod_pow_asm(n as u64, mod_val - 2);
        for x in a.iter_mut() {
            *x = mod_mul(*x, inv_n);
        }
    }
}

pub fn biguint_mul_ntt_mont_asm(a: &BigUInt, b: &BigUInt) -> BigUInt {
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
        fa[i] = montgomery_mul(fa[i], fb[i]);
    }

    ntt_inplace(&mut fa, NTT_MOD, NTT_ROOT, true);

    let mut result = Vec::with_capacity(result_len);
    for i in 0..result_len {
        result.push(fa[i]);
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
        let c = biguint_mul_ntt_mont_asm(&a, &b);
        let expected = BigUInt::from_u64(12345 * 67890);
        assert_eq!(c, expected);
    }
}

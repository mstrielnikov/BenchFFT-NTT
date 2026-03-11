use crate::bigint::{next_power_of_two, BigUInt};

const M61_MOD: u64 = 2305843009213693951;
const M61_P: u64 = 61;

#[inline]
fn m61_reduce(x: u128) -> u64 {
    let t = ((x >> M61_P) as u64).wrapping_add((x & M61_MOD as u128) as u64);
    let t = t.wrapping_add(t >> M61_P);
    (t & M61_MOD) as u64
}

#[inline]
fn m61_add(a: u64, b: u64) -> u64 {
    let x = a.wrapping_add(b);
    let x = x.wrapping_add(x >> M61_P);
    (x & M61_MOD) as u64
}

#[inline]
fn m61_sub(a: u64, b: u64) -> u64 {
    let x = a.wrapping_add(M61_MOD).wrapping_sub(b);
    let x = x.wrapping_add(x >> M61_P);
    (x & M61_MOD) as u64
}

fn m61_mod_pow(mut a: u64, mut e: u64) -> u64 {
    let mut res = 1;
    while e > 0 {
        if e & 1 == 1 {
            res = m61_reduce((res as u128) * (a as u128));
        }
        a = m61_reduce((a as u128) * (a as u128));
        e >>= 1;
    }
    res
}

#[inline]
fn m61_mul(a: u64, b: u64) -> u64 {
    m61_reduce((a as u128) * (b as u128))
}

fn ntt_inplace_m61(a: &mut [u64], inverse: bool) {
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
        let wlen = m61_mod_pow(3, (M61_MOD - 1) / len as u64);
        let wlen = if inverse {
            m61_mod_pow(wlen, M61_MOD - 2)
        } else {
            wlen
        };

        for i in (0..n).step_by(len) {
            let mut w = 1u64;
            for j in 0..len / 2 {
                let u = a[i + j];
                let v = m61_mul(a[i + j + len / 2], w);
                a[i + j] = m61_add(u, v);
                a[i + j + len / 2] = m61_sub(u, v);
                w = m61_mul(w, wlen);
            }
        }
        len <<= 1;
    }

    if inverse {
        let inv_n = m61_mod_pow(n as u64, M61_MOD - 2);
        for x in a.iter_mut() {
            *x = m61_mul(*x, inv_n);
        }
    }
}

pub fn biguint_mul_ntt_mersenne(a: &BigUInt, b: &BigUInt) -> BigUInt {
    if a.is_zero() || b.is_zero() {
        return BigUInt::new();
    }

    let n = next_power_of_two(a.words().len() + b.words().len() - 1);

    let mut fa: Vec<u64> = vec![0; n];
    let mut fb: Vec<u64> = vec![0; n];

    for i in 0..a.words().len() {
        fa[i] = a.words()[i] % M61_MOD;
    }
    for i in 0..b.words().len() {
        fb[i] = b.words()[i] % M61_MOD;
    }

    ntt_inplace_m61(&mut fa, false);
    ntt_inplace_m61(&mut fb, false);

    for i in 0..n {
        fa[i] = m61_mul(fa[i], fb[i]);
    }

    ntt_inplace_m61(&mut fa, true);

    BigUInt::from_slice(&fa)
}

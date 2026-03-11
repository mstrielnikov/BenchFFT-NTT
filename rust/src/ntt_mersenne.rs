use crate::bigint::{next_power_of_two, BigUInt};

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

const M61_MOD: u64 = 2305843009213693951;
const M61_P: u64 = 61;

#[cfg(target_arch = "x86_64")]
#[inline]
unsafe fn m61_reduce_avx(x: __m256i) -> __m256i {
    let mod_mask = _mm256_set1_epi64x(M61_MOD as i64);

    let t = _mm256_srli_epi64(x, M61_P as i32);
    let t = _mm256_add_epi64(t, _mm256_and_si256(x, mod_mask));
    let t = _mm256_add_epi64(t, _mm256_srli_epi64(t, M61_P as i32));
    _mm256_and_si256(t, mod_mask)
}

#[cfg(target_arch = "x86_64")]
#[inline]
unsafe fn m61_add_avx(a: __m256i, b: __m256i) -> __m256i {
    let mod_mask = _mm256_set1_epi64x(M61_MOD as i64);
    let x = _mm256_add_epi64(a, b);
    let x = _mm256_add_epi64(x, _mm256_srli_epi64(x, M61_P as i32));
    _mm256_and_si256(x, mod_mask)
}

#[cfg(target_arch = "x86_64")]
#[inline]
unsafe fn m61_sub_avx(a: __m256i, b: __m256i) -> __m256i {
    let mod_mask = _mm256_set1_epi64x(M61_MOD as i64);
    let x = _mm256_add_epi64(a, mod_mask);
    let x = _mm256_sub_epi64(x, b);
    let x = _mm256_add_epi64(x, _mm256_srli_epi64(x, M61_P as i32));
    _mm256_and_si256(x, mod_mask)
}

#[cfg(target_arch = "x86_64")]
#[inline]
unsafe fn m61_mul_avx(a: __m256i, b: __m256i) -> __m256i {
    let mod_mask = _mm256_set1_epi64x(M61_MOD as i64);

    let a_lo = _mm256_and_si256(a, mod_mask);
    let a_hi = _mm256_srli_epi64(a, M61_P as i32);
    let b_lo = _mm256_and_si256(b, mod_mask);
    let b_hi = _mm256_srli_epi64(b, M61_P as i32);

    let lo_lo = _mm256_mul_epu32(a_lo, b_lo);
    let lo_hi = _mm256_mul_epu32(a_lo, b_hi);
    let hi_lo = _mm256_mul_epu32(a_hi, b_lo);
    let hi_hi = _mm256_mul_epu32(a_hi, b_hi);

    let t1 = _mm256_or_si256(lo_lo, _mm256_slli_epi64(lo_hi, 32));
    let t2 = _mm256_or_si256(hi_lo, _mm256_slli_epi64(hi_hi, 32));
    let prod = _mm256_add_epi64(t1, _mm256_slli_epi64(t2, M61_P as i32));

    m61_reduce_avx(prod)
}

#[cfg(target_arch = "x86_64")]
unsafe fn ntt_inplace_m61_avx(a: &mut [u64], inverse: bool) {
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

    let mut roots: Vec<u64> = vec![0; n];
    roots[0] = 1;
    for len in (2..=n).step_by(2) {
        let wlen = m61_mod_pow(3, (M61_MOD - 1) / len as u64);
        let wlen = if inverse {
            m61_mod_pow(wlen, M61_MOD - 2)
        } else {
            wlen
        };
        for i in len / 2..len {
            roots[i] = m61_mul(roots[i - len / 2], wlen);
        }
    }

    let mut len = 2;
    while len <= n {
        let avx_len = len / 2;
        for i in (0..n).step_by(len) {
            let mut j = 0;
            while j + 8 <= avx_len {
                let idx = i + j;
                let idx2 = i + j + avx_len;

                let u0 = _mm256_loadu_si256(&a[idx] as *const u64 as *const __m256i);
                let u1 = _mm256_loadu_si256(&a[idx + 4] as *const u64 as *const __m256i);
                let v0 = m61_mul_avx(
                    _mm256_loadu_si256(&a[idx2] as *const u64 as *const __m256i),
                    _mm256_loadu_si256(&roots[j] as *const u64 as *const __m256i),
                );
                let v1 = m61_mul_avx(
                    _mm256_loadu_si256(&a[idx2 + 4] as *const u64 as *const __m256i),
                    _mm256_loadu_si256(&roots[j + 4] as *const u64 as *const __m256i),
                );

                let r0 = m61_add_avx(u0, v0);
                let r1 = m61_add_avx(u1, v1);
                let s0 = m61_sub_avx(u0, v0);
                let s1 = m61_sub_avx(u1, v1);

                _mm256_storeu_si256(&mut a[idx] as *mut u64 as *mut __m256i, r0);
                _mm256_storeu_si256(&mut a[idx + 4] as *mut u64 as *mut __m256i, r1);
                _mm256_storeu_si256(&mut a[idx2] as *mut u64 as *mut __m256i, s0);
                _mm256_storeu_si256(&mut a[idx2 + 4] as *mut u64 as *mut __m256i, s1);

                j += 8;
            }
            while j + 4 <= avx_len {
                let idx = i + j;
                let idx2 = i + j + avx_len;

                let u0 = _mm256_loadu_si256(&a[idx] as *const u64 as *const __m256i);
                let v0 = m61_mul_avx(
                    _mm256_loadu_si256(&a[idx2] as *const u64 as *const __m256i),
                    _mm256_loadu_si256(&roots[j] as *const u64 as *const __m256i),
                );

                let r0 = m61_add_avx(u0, v0);
                let s0 = m61_sub_avx(u0, v0);

                _mm256_storeu_si256(&mut a[idx] as *mut u64 as *mut __m256i, r0);
                _mm256_storeu_si256(&mut a[idx2] as *mut u64 as *mut __m256i, s0);

                j += 4;
            }
            while j < avx_len {
                let idx = i + j;
                let idx2 = i + j + avx_len;
                let u = a[idx];
                let v = m61_mul(a[idx2], roots[j]);
                a[idx] = m61_add(u, v);
                a[idx2] = m61_sub(u, v);
                j += 1;
            }
        }
        len <<= 1;
    }

    if inverse {
        let inv_n = m61_mod_pow(n as u64, M61_MOD - 2);
        let inv_n_vec = _mm256_set1_epi64x(inv_n as i64);
        let mut i = 0;
        while i + 4 <= n {
            let val = _mm256_loadu_si256(&a[i] as *const u64 as *const __m256i);
            let val = m61_mul_avx(val, inv_n_vec);
            _mm256_storeu_si256(&mut a[i] as *mut u64 as *mut __m256i, val);
            i += 4;
        }
        while i < n {
            a[i] = m61_mul(a[i], inv_n);
            i += 1;
        }
    }
}

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

#[cfg(target_arch = "x86_64")]
pub fn biguint_mul_ntt_mersenne_avx(a: &BigUInt, b: &BigUInt) -> BigUInt {
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

    unsafe {
        ntt_inplace_m61_avx(&mut fa, false);
        ntt_inplace_m61_avx(&mut fb, false);
    }

    #[cfg(target_arch = "x86_64")]
    unsafe {
        let mut i = 0;
        while i + 4 <= n {
            let a_vec = _mm256_loadu_si256(&fa[i] as *const u64 as *const __m256i);
            let b_vec = _mm256_loadu_si256(&fb[i] as *const u64 as *const __m256i);
            let prod = m61_mul_avx(a_vec, b_vec);
            _mm256_storeu_si256(&mut fa[i] as *mut u64 as *mut __m256i, prod);
            i += 4;
        }
        while i < n {
            fa[i] = m61_mul(fa[i], fb[i]);
            i += 1;
        }
    }

    #[cfg(not(target_arch = "x86_64"))]
    for i in 0..n {
        fa[i] = m61_mul(fa[i], fb[i]);
    }

    unsafe {
        ntt_inplace_m61_avx(&mut fa, true);
    }

    BigUInt::from_slice(&fa)
}

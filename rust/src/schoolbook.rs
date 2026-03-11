use crate::bigint::BigUInt;

pub fn biguint_mul_schoolbook(a: &BigUInt, b: &BigUInt) -> BigUInt {
    if a.is_zero() || b.is_zero() {
        return BigUInt::new();
    }

    let result_len = a.words().len() + b.words().len();
    let mut result = vec![0u64; result_len];

    for i in 0..a.words().len() {
        let mut carry: u128 = 0;
        for j in 0..b.words().len() {
            let prod =
                (a.words()[i] as u128) * (b.words()[j] as u128) + (result[i + j] as u128) + carry;
            result[i + j] = prod as u64;
            carry = prod >> 64;
        }
        if carry != 0 {
            result[i + b.words().len()] += carry as u64;
        }
    }

    BigUInt::from_slice(&result)
}

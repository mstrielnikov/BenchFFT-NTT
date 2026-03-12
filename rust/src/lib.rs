pub mod bigint;
pub mod fft;
pub mod fft_avx;
pub mod fft_mersenne;
pub mod ntt_mersenne;
pub mod ntt_mersenne_alt;
pub mod ntt_mont;
pub mod ntt_mont_asm;
pub mod schoolbook;

pub use bigint::{next_power_of_two, BigUInt};
pub use fft::biguint_mul_fft_split;
pub use fft_avx::biguint_mul_fft_split_avx;
pub use fft_mersenne::biguint_mul_fft_mersenne;
pub use ntt_mersenne::biguint_mul_ntt_mersenne;
pub use ntt_mersenne_alt::biguint_mul_ntt_mersenne_alt;
pub use ntt_mont::biguint_mul_ntt_mont;
pub use ntt_mont_asm::biguint_mul_ntt_mont_asm;
pub use schoolbook::biguint_mul_schoolbook;

#[cfg(target_arch = "x86_64")]
pub use ntt_mersenne::biguint_mul_ntt_mersenne_avx;

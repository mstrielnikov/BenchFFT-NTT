use benchfft::{
    biguint_mul_fft_split, biguint_mul_ntt_mersenne, biguint_mul_ntt_mont, biguint_mul_schoolbook,
    BigUInt,
};

#[cfg(target_arch = "x86_64")]
use benchfft::{biguint_mul_fft_split_avx, biguint_mul_ntt_mersenne_avx};

use std::time::Instant;

#[allow(dead_code)]
fn run_benchmark(_name: &str, size: usize, iterations: usize) {
    let mut seed: u64 = 42;
    let words1: Vec<u64> = (0..size)
        .map(|_| {
            seed = seed.wrapping_mul(1103515245).wrapping_add(12345);
            seed
        })
        .collect();
    let words2: Vec<u64> = (0..size)
        .map(|_| {
            seed = seed.wrapping_mul(1103515245).wrapping_add(12345);
            seed
        })
        .collect();

    let a = BigUInt::from_slice(&words1);
    let b = BigUInt::from_slice(&words2);

    let start = Instant::now();
    for _ in 0..iterations {
        let _c = biguint_mul_fft_split(&a, &b);
    }
    let elapsed = start.elapsed().as_secs_f64() * 1000.0;
    println!("FFT {:4}:   {:.4} ms ({} iters)", size, elapsed, iterations);

    let start = Instant::now();
    for _ in 0..iterations {
        let _c = biguint_mul_ntt_mont(&a, &b);
    }
    let elapsed = start.elapsed().as_secs_f64() * 1000.0;
    println!("NTT {:4}:   {:.4} ms ({} iters)", size, elapsed, iterations);
}

#[cfg(target_arch = "x86_64")]
fn run_benchmark_avx(_name: &str, size: usize, iterations: usize) {
    let mut seed: u64 = 42;
    let words1: Vec<u64> = (0..size)
        .map(|_| {
            seed = seed.wrapping_mul(1103515245).wrapping_add(12345);
            seed
        })
        .collect();
    let words2: Vec<u64> = (0..size)
        .map(|_| {
            seed = seed.wrapping_mul(1103515245).wrapping_add(12345);
            seed
        })
        .collect();

    let a = BigUInt::from_slice(&words1);
    let b = BigUInt::from_slice(&words2);

    let start = Instant::now();
    for _ in 0..iterations {
        let _c = biguint_mul_fft_split_avx(&a, &b);
    }
    let elapsed = start.elapsed().as_secs_f64() * 1000.0;
    println!("FFT {:4}:   {:.4} ms ({} iters)", size, elapsed, iterations);

    let start = Instant::now();
    for _ in 0..iterations {
        let _c = biguint_mul_ntt_mont(&a, &b);
    }
    let elapsed = start.elapsed().as_secs_f64() * 1000.0;
    println!("NTT {:4}:   {:.4} ms ({} iters)", size, elapsed, iterations);
}

fn run_benchmark_mersenne(_name: &str, size: usize, iterations: usize) {
    let mut seed: u64 = 42;
    let words1: Vec<u64> = (0..size)
        .map(|_| {
            seed = seed.wrapping_mul(1103515245).wrapping_add(12345);
            seed
        })
        .collect();
    let words2: Vec<u64> = (0..size)
        .map(|_| {
            seed = seed.wrapping_mul(1103515245).wrapping_add(12345);
            seed
        })
        .collect();

    let a = BigUInt::from_slice(&words1);
    let b = BigUInt::from_slice(&words2);

    let start = Instant::now();
    for _ in 0..iterations {
        let _c = biguint_mul_ntt_mersenne(&a, &b);
    }
    let elapsed = start.elapsed().as_secs_f64() * 1000.0;
    println!(
        "Mersenne NTT {:4}:   {:.4} ms ({} iters)",
        size, elapsed, iterations
    );
}

#[cfg(target_arch = "x86_64")]
fn run_benchmark_mersenne_avx(_name: &str, size: usize, iterations: usize) {
    let mut seed: u64 = 42;
    let words1: Vec<u64> = (0..size)
        .map(|_| {
            seed = seed.wrapping_mul(1103515245).wrapping_add(12345);
            seed
        })
        .collect();
    let words2: Vec<u64> = (0..size)
        .map(|_| {
            seed = seed.wrapping_mul(1103515245).wrapping_add(12345);
            seed
        })
        .collect();

    let a = BigUInt::from_slice(&words1);
    let b = BigUInt::from_slice(&words2);

    let start = Instant::now();
    for _ in 0..iterations {
        let _c = biguint_mul_ntt_mersenne_avx(&a, &b);
    }
    let elapsed = start.elapsed().as_secs_f64() * 1000.0;
    println!(
        "Mersenne NTT AVX {:4}:   {:.4} ms ({} iters)",
        size, elapsed, iterations
    );
}

fn run_benchmark_schoolbook(_name: &str, size: usize, iterations: usize) {
    let mut seed: u64 = 42;
    let words1: Vec<u64> = (0..size)
        .map(|_| {
            seed = seed.wrapping_mul(1103515245).wrapping_add(12345);
            seed
        })
        .collect();
    let words2: Vec<u64> = (0..size)
        .map(|_| {
            seed = seed.wrapping_mul(1103515245).wrapping_add(12345);
            seed
        })
        .collect();

    let a = BigUInt::from_slice(&words1);
    let b = BigUInt::from_slice(&words2);

    let start = Instant::now();
    for _ in 0..iterations {
        let _c = biguint_mul_schoolbook(&a, &b);
    }
    let elapsed = start.elapsed().as_secs_f64() * 1000.0;
    println!(
        "Schoolbook {:4}:   {:.4} ms ({} iters)",
        size, elapsed, iterations
    );
}

fn main() {
    #[cfg(target_arch = "x86_64")]
    {
        println!("========== RUST BENCHMARKS (Native CPU + AVX Intrinsics) ==========\n");
        println!("ML-KEM / ML-DSA Vector Sizes");
        println!("=============================\n");

        println!("ML-KEM-512 (256 words):");
        run_benchmark_avx("FFT", 256, 100);
        println!();

        println!("ML-KEM-768 (512 words):");
        run_benchmark_avx("FFT", 512, 50);
        println!();

        println!("ML-KEM-1024 (1024 words):");
        run_benchmark_avx("FFT", 1024, 20);
        println!();

        println!("ML-DSA (2048 words):");
        run_benchmark_avx("FFT", 2048, 10);
        println!();

        println!("Extended (3072 words):");
        run_benchmark_avx("FFT", 3072, 5);
        println!();

        println!("Extended (4096 words):");
        run_benchmark_avx("FFT", 4096, 3);
        println!();

        println!("================================\n");
    }

    println!("========== RUST MERSENNE NTT BENCHMARKS ==========\n");
    println!("ML-KEM / ML-DSA Vector Sizes\n");

    println!("ML-KEM-512 (256 words):");
    run_benchmark_mersenne("Mersenne", 256, 100);
    #[cfg(target_arch = "x86_64")]
    run_benchmark_mersenne_avx("Mersenne", 256, 100);
    println!();

    println!("ML-KEM-768 (512 words):");
    run_benchmark_mersenne("Mersenne", 512, 50);
    #[cfg(target_arch = "x86_64")]
    run_benchmark_mersenne_avx("Mersenne", 512, 50);
    println!();

    println!("ML-KEM-1024 (1024 words):");
    run_benchmark_mersenne("Mersenne", 1024, 20);
    #[cfg(target_arch = "x86_64")]
    run_benchmark_mersenne_avx("Mersenne", 1024, 20);
    println!();

    println!("ML-DSA (2048 words):");
    run_benchmark_mersenne("Mersenne", 2048, 10);
    #[cfg(target_arch = "x86_64")]
    run_benchmark_mersenne_avx("Mersenne", 2048, 10);
    println!();

    println!("Extended (3072 words):");
    run_benchmark_mersenne("Mersenne", 3072, 5);
    #[cfg(target_arch = "x86_64")]
    run_benchmark_mersenne_avx("Mersenne", 3072, 5);
    println!();

    println!("Extended (4096 words):");
    run_benchmark_mersenne("Mersenne", 4096, 3);
    #[cfg(target_arch = "x86_64")]
    run_benchmark_mersenne_avx("Mersenne", 4096, 3);
    println!();

    println!("================================\n");

    println!("========== RUST SCHOOLBOOK BENCHMARKS ==========\n");
    println!("ML-KEM / ML-DSA Vector Sizes\n");

    println!("ML-KEM-512 (256 words):");
    run_benchmark_schoolbook("Schoolbook", 256, 100);
    println!();

    println!("ML-KEM-768 (512 words):");
    run_benchmark_schoolbook("Schoolbook", 512, 50);
    println!();

    println!("ML-KEM-1024 (1024 words):");
    run_benchmark_schoolbook("Schoolbook", 1024, 20);
    println!();

    println!("ML-DSA (2048 words):");
    run_benchmark_schoolbook("Schoolbook", 2048, 10);
    println!();

    println!("Extended (3072 words):");
    run_benchmark_schoolbook("Schoolbook", 3072, 5);
    println!();

    println!("Extended (4096 words):");
    run_benchmark_schoolbook("Schoolbook", 4096, 3);
    println!();

    println!("================================");
}

use benchfft::{biguint_mul_fft_split, biguint_mul_ntt_mont, BigUInt};
use std::time::Instant;

fn run_benchmark(name: &str, size: usize, iterations: usize) {
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

    // FFT Benchmark
    let start = Instant::now();
    for _ in 0..iterations {
        let _c = biguint_mul_fft_split(&a, &b);
    }
    let elapsed = start.elapsed().as_secs_f64() * 1000.0;
    println!("FFT {:4}:   {:.4} ms ({} iters)", size, elapsed, iterations);

    // NTT Benchmark
    let start = Instant::now();
    for _ in 0..iterations {
        let _c = biguint_mul_ntt_mont(&a, &b);
    }
    let elapsed = start.elapsed().as_secs_f64() * 1000.0;
    println!("NTT {:4}:   {:.4} ms ({} iters)", size, elapsed, iterations);
}

fn main() {
    println!("========== RUST BENCHMARKS ==========\n");
    println!("ML-KEM / ML-DSA Vector Sizes");
    println!("=============================\n");

    // ML-KEM sizes (k = 2,3,4)
    println!("ML-KEM-512 (256 words):");
    run_benchmark("FFT", 256, 100);
    println!();

    println!("ML-KEM-768 (512 words):");
    run_benchmark("FFT", 512, 50);
    println!();

    println!("ML-KEM-1024 (1024 words):");
    run_benchmark("FFT", 1024, 20);
    println!();

    // ML-DSA sizes
    println!("ML-DSA (2048 words):");
    run_benchmark("FFT", 2048, 10);
    println!();

    // Extended sizes
    println!("Extended (3072 words):");
    run_benchmark("FFT", 3072, 5);
    println!();

    println!("Extended (4096 words):");
    run_benchmark("FFT", 4096, 3);
    println!();

    println!("================================");
}

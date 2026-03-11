# BenchFFT-NTT

Big integer multiplication using FFT and NTT convolution algorithms. Mirror implementations in C and Rust for benchmarking.

## Project Structure

```
BenchFFT-NTT/
├── c/                      # C implementation
│   ├── Makefile
│   ├── include/bigint.h   # Header
│   ├── src/
│   │   ├── bigint.c       # Core BigUInt
│   │   ├── fft_split.c    # FFT implementation
│   │   ├── fft_split_avx.c # FFT with AVX intrinsics
│   │   ├── ntt_mont.c     # NTT Montgomery (998244353)
│   │   ├── ntt_mersenne.c # NTT Mersenne (2^61-1)
│   │   └── ntt_mersenne_avx.c # NTT Mersenne + AVX
│   └── test/
│       ├── test.c          # Tests
│       └── bench.c         # Benchmarks
├── rust/                   # Rust implementation
│   ├── Cargo.toml
│   └── src/
│       ├── lib.rs         # Module re-exports
│       ├── bigint.rs      # BigUInt type
│       ├── fft.rs         # Basic FFT
│       ├── fft_avx.rs     # AVX-optimized FFT
│       ├── ntt_mont.rs    # Montgomery NTT
│       ├── ntt_mersenne.rs # Mersenne NTT + AVX
│       └── main.rs        # Benchmark runner
└── lean/                  # Lean4 formal verification (placeholder)
    ├── Formal.lean
    └── README.md
```

## Algorithm Implementations

### C Modules

| File | Description |
|------|-------------|
| `bigint.c` | Core BigUInt operations |
| `fft_split.c` | Basic Cooley-Tukey FFT |
| `fft_split_avx.c` | AVX2 vectorized FFT |
| `ntt_mont.c` | Montgomery NTT (mod 998244353) |
| `ntt_mersenne.c` | Mersenne NTT (mod 2^61-1) |
| `ntt_mersenne_avx.c` | Mersenne NTT with AVX |

### Rust Modules

| File | Description |
|------|-------------|
| `bigint.rs` | BigUInt type definition |
| `fft.rs` | Basic Cooley-Tukey FFT |
| `fft_avx.rs` | AVX2 vectorized FFT |
| `ntt_mont.rs` | Montgomery NTT (mod 998244353) |
| `ntt_mersenne.rs` | Mersenne NTT + AVX (mod 2^61-1) |

## Algorithms

### Fast Fourier Transform (FFT)

The Cooley-Tukey FFT algorithm computes the Discrete Fourier Transform (DFT) in O(n log n) time instead of O(n²). For polynomial multiplication:

1. Convert coefficients to complex numbers
2. Apply FFT to both polynomials
3. Pointwise multiply the transforms
4. Apply inverse FFT to get result coefficients

### Number Theoretic Transform (NTT)

NTT is the FFT over a finite field instead of complex numbers. Two moduli are tested:

#### Montgomery (998244353)
- **Modulus**: 998244353 (2^23 * 7 * 17 + 1)
- **Primitive Root**: 3
- Standard NTT using modular arithmetic with `%` operator

#### Mersenne (2^61 - 1)
- **Modulus**: 2305843009213693951 (2^61 - 1)
- **Reduction**: Uses shift-and-add instead of division
- **Formula**: `x mod M61 = (x & M61) + (x >> 61)`

Advantages of Mersenne:
- No expensive division operations
- Uses simple bit shifts
- ~3-5x faster than Montgomery

## Building

```bash
# Build C library
make build                    # Default (auto-vectorized)
make build-scalar             # No vectorization
make build-avx               # AVX intrinsics
make build-mersenne          # Mersenne NTT
make build-mersenne-avx      # Mersenne + AVX

# Run C benchmarks
make bench
make bench-scalar
make bench-avx
make bench-mersenne
make bench-mersenne-avx

# Rust benchmarks
cd rust && cargo build --release
cargo run --release

# Clean
make clean
```

## Benchmark Results

All times in milliseconds. Lower is better.

### Test Environment
- **CPU**: Native x86_64 with AVX2/FMA
- **Compiler**: GCC (C), Rust 1.75+ (Rust)
- **Optimization**: -O3 with target-cpu=native

---

## C Implementation Benchmarks

### FFT (Complex Number)

| Size | Auto-Vectorized | AVX Intrinsics | Scalar |
|------|-----------------|----------------|--------|
| 256  | 4.96 ms | 5.05 ms | 8.77 ms |
| 512  | 2.49 ms | 1.65 ms | 5.24 ms |
| 1024 | 2.09 ms | 1.37 ms | 2.14 ms |
| 2048 | 2.47 ms | 1.56 ms | 2.44 ms |
| 3072 | 3.63 ms | 2.77 ms | 3.74 ms |
| 4096 | 1.76 ms | 1.31 ms | 1.80 ms |

### NTT Montgomery (mod 998244353)

| Size | Auto-Vectorized | AVX Intrinsics | Scalar |
|------|-----------------|----------------|--------|
| 256  | 2.22 ms | 2.22 ms | 2.37 ms |
| 512  | 7.57 ms | 7.46 ms | 7.59 ms |
| 1024 | 6.24 ms | 6.29 ms | 6.31 ms |
| 2048 | 6.63 ms | 6.57 ms | 6.79 ms |
| 3072 | 6.99 ms | 7.00 ms | 7.11 ms |
| 4096 | 4.17 ms | 4.22 ms | 4.20 ms |

### NTT Mersenne (mod 2^61-1)

| Size | Scalar | AVX2 Intrinsics |
|------|--------|-----------------|
| 256  | 4.99 ms | 2.22 ms |
| 512  | 3.55 ms | 1.94 ms |
| 1024 | 1.53 ms | 1.61 ms |
| 2048 | 1.56 ms | 1.63 ms |
| 3072 | 1.64 ms | 1.85 ms |
| 4096 | 1.00 ms | 1.01 ms |

---

## Rust Implementation Benchmarks

### FFT (Complex Number)

| Size | Scalar | AVX Intrinsics |
|------|--------|----------------|
| 256  | 20.87 ms | 16.76 ms |
| 512  | 14.08 ms | 13.87 ms |
| 1024 | 12.44 ms | 12.31 ms |
| 2048 | 13.94 ms | 13.91 ms |
| 3072 | 15.00 ms | 15.73 ms |
| 4096 | 10.01 ms | 9.39 ms |

### NTT Montgomery (mod 998244353)

| Size | Scalar |
|------|--------|
| 256  | 10.96 ms |
| 512  | 11.50 ms |
| 1024 | 9.90 ms |
| 2048 | 10.63 ms |
| 3072 | 11.75 ms |
| 4096 | 7.05 ms |

### NTT Mersenne (mod 2^61-1)

| Size | Scalar | AVX2 Intrinsics |
|------|--------|------------------|
| 256  | 4.13 ms | 68.73 ms |
| 512  | 4.12 ms | 88.11 ms |
| 1024 | 3.52 ms | 95.25 ms |
| 2048 | 3.67 ms | 140.75 ms |
| 3072 | 3.72 ms | 228.70 ms |
| 4096 | 2.29 ms | 137.00 ms |

### Schoolbook Multiplication

| Size | Time |
|------|------|
| 256  | 6.26 ms |
| 512  | 12.61 ms |
| 1024 | 19.98 ms |
| 2048 | 39.98 ms |
| 3072 | 44.93 ms |
| 4096 | 47.79 ms |

---

## Cross-Platform Comparison

### 4096-word multiplication

| Implementation | FFT | NTT Montgomery | NTT Mersenne | Schoolbook |
|----------------|-----|----------------|--------------|------------|
| **C Auto** | 1.76 ms | 4.17 ms | 1.00 ms | - |
| **C AVX** | 1.31 ms | 4.22 ms | 1.01 ms | - |
| **C Scalar** | 1.80 ms | 4.20 ms | - | - |
| **Rust Scalar** | 10.01 ms | 7.05 ms | 2.29 ms | 47.79 ms |
| **Rust AVX** | 9.39 ms | - | N/A | - |

### 1024-word multiplication

| Implementation | FFT | NTT Montgomery | NTT Mersenne |
|----------------|-----|----------------|--------------|
| **C Auto** | 2.09 ms | 6.24 ms | 1.53 ms |
| **C AVX** | 1.37 ms | 6.29 ms | 1.61 ms |
| **Rust Scalar** | 12.44 ms | 9.90 ms | 3.52 ms |
| **Rust AVX** | 12.31 ms | - | N/A |

---

## Performance Rankings

### Best for Large Multiplication (4096 words)

| Rank | Algorithm | Language | Time |
|------|-----------|----------|------|
| 1 | C Mersenne NTT | C | 1.00 ms |
| 2 | C Mersenne AVX | C | 1.01 ms |
| 3 | C FFT AVX | C | 1.31 ms |
| 4 | C FFT Auto | C | 1.76 ms |
| 5 | Rust Mersenne | Rust | 2.29 ms |
| 6 | C FFT Scalar | C | 1.80 ms |
| 7 | C NTT Montgomery | C | 4.17 ms |
| 8 | Rust FFT AVX | Rust | 9.39 ms |
| 9 | Rust NTT Montgomery | Rust | 7.05 ms |

---

## Key Findings

### 1. Mersenne > Montgomery: 3-5x faster
The shift-add reduction in Mersenne (2^61-1) avoids expensive modulo division required by Montgomery (998244353).

### 2. C > Rust: 2-5x faster
Even with equivalent algorithms, C outperforms Rust significantly:
- C Mersenne: 1.00 ms vs Rust: 2.29 ms (2.3x faster)
- C FFT: 1.31 ms vs Rust: 9.39 ms (7.2x faster)

### 3. AVX Intrinsics vs Auto-Vectorization
- **FFT**: AVX intrinsics provide 20-35% improvement
- **Montgomery NTT**: Minimal difference (auto-vectorization handles well)
- **Mersenne NTT**: AVX version actually slower than scalar (M61 arithmetic doesn't map well to SIMD)

### 4. SIMD Unfriendly: Mersenne AVX
The Mersenne AVX implementations are significantly slower than scalar because:
- 61-bit arithmetic doesn't fit neatly into 64-bit SIMD lanes
- Requires complex emulation with 32-bit multiplies
- Native 128-bit integer (u128) in scalar code is more efficient

### 5. FFT vs NTT
- **FFT**: Good accuracy with floating-point, ~1.3-1.8 ms for 4096
- **Mersenne NTT**: Fastest at 1.0 ms, exact arithmetic
- **Montgomery NTT**: Slowest at 4.2 ms, requires division

---

## Compilation Flags

| Config | Flags |
|--------|-------|
| C Scalar | `-O3 -std=c11 -fno-tree-vectorize` |
| C Auto-vectorized | `-O3 -std=c11 -mavx2 -mfma -funroll-loops` |
| C AVX Intrinsics | `-O3 -std=c11 -mavx2 -mfma -funroll-loops -DBUILD_AVX` |
| Rust | `--release -C target-cpu=native` |

---

## Testing

- **C**: Tests covering sizes 256, 512, 1024, 2048, 3072, 4096
- **Rust**: Unit tests for small and large multiplications

## Formal Verification

Lean4 formal verification setup in `lean/` directory. Uses FFI + property-based testing approach to verify C implementations against reference implementations.

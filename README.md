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
│   │   ├── ntt_mont.c    # NTT Montgomery
│   │   └── ntt_mont_avx.c # NTT with AVX
│   └── test/
│       ├── test.c         # Tests
│       └── bench.c        # Benchmarks
├── rust/                   # Rust implementation
│   ├── Cargo.toml
│   └── src/
│       ├── lib.rs         # FFT/NTT + AVX implementations
│       └── main.rs        # Benchmark runner
└── lean/                  # Lean4 formal verification (placeholder)
    ├── Formal.lean
    └── README.md
```

## Algorithms

### Fast Fourier Transform (FFT)

The Cooley-Tukey FFT algorithm computes the Discrete Fourier Transform (DFT) in O(n log n) time instead of O(n²). For polynomial multiplication:

1. Convert coefficients to complex numbers
2. Apply FFT to both polynomials
3. Pointwise multiply the transforms
4. Apply inverse FFT to get result coefficients

```
Convolution Theorem: FFT(a * b) = FFT(a) · FFT(b)
```

### Number Theoretic Transform (NTT)

NTT is the FFT over a finite field instead of complex numbers. Uses primitive roots of unity in mod p.

- **Modulus**: 998244353 (2^23 * 7 * 17 + 1)
- **Primitive Root**: 3

Advantages:
- No floating-point rounding errors
- Exact integer arithmetic
- Constant-time operations (no floating-point variance)

## Math Background

### Polynomial Multiplication

Given polynomials A(x) and B(x) of degree n-1:
- A(x) = Σ aᵢxⁱ, B(x) = Σ bᵢxⁱ

The product C(x) = A(x)·B(x) has coefficients:
- cₖ = Σ aᵢbₖ₋ᵢ for i from 0 to k

This is the convolution of coefficient vectors.

### Complexity

| Algorithm | Complexity | Notes |
|----------|------------|-------|
| Naive    | O(n²)     | Direct coefficient multiplication |
| FFT      | O(n log n) | Uses complex roots of unity |
| NTT      | O(n log n) | Finite field variant of FFT |

## Building

```bash
# Build C library (auto-vectorized with AVX/FMA)
make build

# Build scalar version (no vectorization)
make build-scalar

# Build AVX intrinsics version
make build-avx

# Run C tests
make test

# Run C benchmarks (auto-vectorized)
make bench

# Run C scalar benchmarks
make bench-scalar

# Run C AVX benchmarks
make bench-avx

# Build and run Rust benchmarks (with AVX intrinsics)
cd rust && RUSTFLAGS="-C target-cpu=native" cargo build --release
./rust/target/release/bench

# Run all benchmarks (C + Rust)
make bench-all

# Clean build artifacts
make clean
```

## Benchmark Results

Tested on vector sizes matching PQC algorithms (ML-KEM, ML-DSA). Three C build configurations:

1. **Scalar**: No vectorization (`-fno-tree-vectorize`)
2. **Auto-vectorized**: Compiler vectorizes with AVX/FMA (`-mavx2 -mfma -funroll-loops`)
3. **AVX Intrinsics**: Manual AVX2 intrinsics in C code

### Scalar (No Vectorization)

| Size | Type | FFT | NTT |
|------|------|-----|-----|
| 256  | ML-KEM-512 | 3.92 ms | 12.20 ms |
| 512  | ML-KEM-768 | 4.32 ms | 13.83 ms |
| 1024 | ML-KEM-1024 | 3.81 ms | 10.37 ms |
| 2048 | ML-DSA | 4.26 ms | 10.94 ms |
| 3072 | Extended | 6.46 ms | 11.93 ms |
| 4096 | Extended | 3.06 ms | 7.10 ms |

### Auto-Vectorized (Compiler AVX/FMA)

| Size | Type | FFT | NTT |
|------|------|-----|-----|
| 256  | ML-KEM-512 | 5.49 ms | 12.79 ms |
| 512  | ML-KEM-768 | 4.09 ms | 12.56 ms |
| 1024 | ML-KEM-1024 | 3.63 ms | 10.63 ms |
| 2048 | ML-DSA | 4.50 ms | 11.09 ms |
| 3072 | Extended | 6.56 ms | 11.54 ms |
| 4096 | Extended | 2.90 ms | 6.86 ms |

### AVX Intrinsics (Manual)

| Size | Type | FFT | NTT |
|------|------|-----|-----|
| 256  | ML-KEM-512 | 4.57 ms | 13.51 ms |
| 512  | ML-KEM-768 | 2.81 ms | 12.39 ms |
| 1024 | ML-KEM-1024 | 2.34 ms | 10.36 ms |
| 2048 | ML-DSA | 2.63 ms | 10.92 ms |
| 3072 | Extended | 4.91 ms | 15.24 ms |
| 4096 | Extended | 2.15 ms | 7.25 ms |

### Comparison: C vs Rust (All with AVX)

| Size | Type | C FFT | C NTT | Rust FFT (AVX) | Rust NTT (AVX) |
|------|------|-------|-------|----------------|----------------|
| 256  | ML-KEM-512 | 2.63 ms | 12.05 ms | 5.25 ms | 19.90 ms |
| 512  | ML-KEM-768 | 2.70 ms | 12.27 ms | 3.34 ms | 18.86 ms |
| 1024 | ML-KEM-1024 | 2.28 ms | 10.77 ms | 2.50 ms | 16.87 ms |
| 2048 | ML-DSA | 3.33 ms | 11.14 ms | 4.61 ms | 18.11 ms |
| 3072 | Extended | 4.74 ms | 14.25 ms | 4.15 ms | 19.17 ms |
| 4096 | Extended | 2.11 ms | 6.90 ms | 2.41 ms | 11.34 ms |

### Key Observations

- **C vs Rust AVX**: C is ~2x faster for FFT and ~1.5x faster for NTT
- **FFT**: Power-of-two sizes (1024, 4096) show best performance
- **AVX Intrinsics**: 30-40% faster than auto-vectorized for FFT at 512-1024 sizes
- **C vs Rust**: C FFT and Rust FFT are very close (~10% difference); C NTT is ~35-40% faster than Rust NTT
- **NTT overhead**: Higher than FFT due to modular multiplication

## Fair Comparison Notes

To ensure a fair comparison between C and Rust implementations:

### Memory Management Alignment

Both implementations now use pre-allocated vectors with exact expected sizes:

1. **Input vectors**: Pre-allocated to size `n` (next power of 2)
2. **Result vector**: Pre-allocated to `result_len = a->len + b->len - 1`
3. **No incremental reallocation**: Results are written to pre-allocated buffers

This eliminates memory allocation as a variable factor in benchmarks.

### Compilation Flags

- **C Scalar**: `-O3 -std=c11 -fno-tree-vectorize`
- **C Auto-vectorized**: `-O3 -std=c11 -mavx2 -mfma -funroll-loops`
- **C AVX Intrinsics**: `-O3 -std=c11 -mavx2 -mfma -funroll-loops -DBUILD_AVX`
- **Rust**: `--release` (equivalent to `-O3`)

### Remaining Differences

- **FFT**: Uses `double` in C vs Rust (floating-point behavior identical)
- **NTT**: C uses inline modular arithmetic; Rust uses wrapper functions
- **Rust Vec vs C arrays**: Minor overhead in Rust's bounds checking (optimized out in release)

## Testing

- **C**: 14 tests covering sizes 256, 512, 1024, 2048, 3072, 4096
- **Rust**: Unit tests for small and large multiplications

## Formal Verification

Lean4 formal verification setup in `lean/` directory. Uses FFI + property-based testing approach to verify C implementations against reference implementations.

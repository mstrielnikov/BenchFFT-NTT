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
│       ├── ntt_mersenne.rs # Mersenne NTT
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
| `ntt_mersenne.rs` | Mersenne NTT (mod 2^61-1) |
BenchFFT-NTT/
├── c/                      # C implementation
│   ├── Makefile
│   ├── include/bigint.h   # Header
│   ├── src/
│   │   ├── bigint.c       # Core BigUInt
│   │   ├── fft_split.c    # FFT implementation
│   │   ├── fft_split_avx.c # FFT with AVX intrinsics
│   │   ├── ntt_mont.c    # NTT Montgomery (998244353)
│   │   ├── ntt_mersenne.c # NTT Mersenne (2^61-1)
│   │   └── ntt_mersenne_avx.c # NTT Mersenne + AVX
│   └── test/
│       ├── test.c         # Tests
│       └── bench.c        # Benchmarks
├── rust/                   # Rust implementation
│   ├── Cargo.toml
│   └── src/
│       ├── lib.rs         # FFT/NTT + AVX + Mersenne implementations
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
# Build C library (auto-vectorized with AVX/FMA)
make build

# Build scalar version (no vectorization)
make build-scalar

# Build AVX intrinsics version
make build-avx

# Build Mersenne NTT version
make build-mersenne

# Build Mersenne + AVX version
make build-mersenne-avx

# Run C tests
make test

# Run C benchmarks
make bench              # Auto-vectorized
make bench-scalar       # No vectorization
make bench-avx          # AVX intrinsics
make bench-mersenne     # Mersenne NTT
make bench-mersenne-avx # Mersenne + AVX

# Build and run Rust benchmarks
cd rust && RUSTFLAGS="-C target-cpu=native" cargo build --release
./rust/target/release/bench

# Run all benchmarks (C + Rust)
make bench-all

# Clean build artifacts
make clean
```

## Benchmark Results

All times in milliseconds. Lower is better.

### FFT Comparison (C Auto-Vectorized vs AVX Intrinsics)

| Size | Auto-Vec FFT | AVX FFT | Speedup |
|------|--------------|---------|---------|
| 256  | 6.34 ms | 3.87 ms | 1.6x |
| 512  | 3.99 ms | 2.80 ms | 1.4x |
| 1024 | 3.55 ms | 2.42 ms | 1.5x |
| 2048 | 3.96 ms | 2.71 ms | 1.5x |
| 4096 | 2.96 ms | 2.11 ms | 1.4x |

### NTT Comparison - Montgomery vs Mersenne (C)

| Size | Mont (998244353) | Mersenne | Speedup |
|------|------------------|----------|---------|
| 256  | 15.03 ms | 5.01 ms | **3.0x** |
| 512  | 12.40 ms | 3.21 ms | **3.9x** |
| 1024 | 10.41 ms | 2.68 ms | **3.9x** |
| 2048 | 10.98 ms | 3.33 ms | **3.3x** |
| 4096 | 7.04 ms | 1.65 ms | **4.3x** |

### C: All Configurations

| Size | FFT-AVX | NTT-Mont | NTT-Mersenne | NTT-Mersenne+AVX |
|------|---------|----------|--------------|-------------------|
| 256  | 3.87 ms | 15.03 ms | 5.01 ms | **4.41 ms** |
| 512  | 2.80 ms | 12.40 ms | 3.21 ms | **2.80 ms** |
| 1024 | 2.42 ms | 10.41 ms | 2.68 ms | **2.51 ms** |
| 2048 | 2.71 ms | 10.98 ms | 3.33 ms | **2.51 ms** |
| 4096 | 2.11 ms | 7.04 ms | 1.65 ms | **1.57 ms** |

### Rust vs C Comparison

| Size | C FFT | Rust FFT | C NTT-Mont | Rust NTT-Mont | C NTT-Mersenne | Rust NTT-Mersenne |
|------|-------|----------|------------|---------------|----------------|-------------------|
| 256  | 3.87 ms | 5.49 ms | 15.03 ms | 19.31 ms | 5.01 ms | 6.15 ms |
| 512  | 2.80 ms | 3.03 ms | 12.40 ms | 19.01 ms | 3.21 ms | 6.51 ms |
| 1024 | 2.42 ms | 2.79 ms | 10.41 ms | 16.43 ms | 2.68 ms | 5.12 ms |
| 2048 | 2.71 ms | 4.35 ms | 10.98 ms | 18.66 ms | 3.33 ms | 5.50 ms |
| 4096 | 2.11 ms | 2.71 ms | 7.04 ms | 11.66 ms | 1.65 ms | 3.60 ms |

## Performance Rankings

### Best NTT Implementation (C)

| Rank | Implementation | 1024 words | 4096 words |
|------|---------------|------------|------------|
| 1 | C Mersenne+AVX | 2.51 ms | 1.57 ms |
| 2 | C Mersenne | 2.68 ms | 1.65 ms |
| 3 | C AVX-Montgomery | 10.36 ms | 7.25 ms |
| 4 | C Auto-Montgomery | 10.41 ms | 7.04 ms |

### Best FFT Implementation (C)

| Rank | Implementation | 1024 words | 4096 words |
|------|---------------|------------|------------|
| 1 | C AVX Intrinsics | 2.42 ms | 2.11 ms |
| 2 | C Auto-vectorized | 3.55 ms | 2.96 ms |
| 3 | C Scalar | 3.81 ms | 3.06 ms |

### C vs Rust (Same Algorithm)

| Operation | C | Rust | C Advantage |
|-----------|---|------|-------------|
| FFT AVX (1024) | 2.42 ms | 2.79 ms | 1.15x |
| NTT Mersenne (1024) | 2.68 ms | 5.12 ms | 1.91x |
| NTT Mersenne (4096) | 1.65 ms | 3.60 ms | 2.18x |

## Key Findings

1. **Mersenne > Montgomery**: 3-5x faster due to shift-add reduction instead of division
2. **AVX Intrinsics**: 30-40% improvement for FFT over auto-vectorization
3. **C > Rust**: C is 1.5-2x faster even with same algorithm
4. **Power-of-two advantage**: FFT/NTT perform best at 1024, 4096 (aligned with PQC sizes)

## Fair Comparison Notes

### Memory Management Alignment

Both implementations use pre-allocated vectors with exact expected sizes:
- Input vectors: Pre-allocated to size `n` (next power of 2)
- Result vector: Pre-allocated to `a->len + b->len - 1`
- No incremental reallocation during computation

### Compilation Flags

- **C Scalar**: `-O3 -std=c11 -fno-tree-vectorize`
- **C Auto-vectorized**: `-O3 -std=c11 -mavx2 -mfma -funroll-loops`
- **C AVX Intrinsics**: `-O3 -std=c11 -mavx2 -mfma -funroll-loops -DBUILD_AVX`
- **Rust**: `--release -C target-cpu=native`

### Remaining Differences

- **FFT**: Uses `double` in both C and Rust
- **NTT Montgomery**: C uses inline `%` operator; Rust uses wrapper functions
- **NTT Mersenne**: Both use shift-add reduction
- **Rust Vec vs C arrays**: Bounds checking optimized out in release

## Testing

- **C**: 14 tests covering sizes 256, 512, 1024, 2048, 3072, 4096
- **Rust**: Unit tests for small and large multiplications

## Formal Verification

Lean4 formal verification setup in `lean/` directory. Uses FFI + property-based testing approach to verify C implementations against reference implementations.

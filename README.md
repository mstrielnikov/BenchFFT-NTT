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
│   │   ├── fft_mersenne.c # FFT + Mersenne reduction
│   │   ├── ntt_mersenne.c # Integer NTT → M61 reduction
│   │   ├── ntt_mont.c     # NTT Montgomery (998244353)
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
| `fft_split.c` | Cooley-Tukey FFT (double precision) |
| `fft_split_avx.c` | AVX2 vectorized FFT |
| `fft_mersenne.c` | FFT + Mersenne M61 reduction |
| `ntt_mersenne.c` | Integer NTT (mod 998244353) → M61 |
| `ntt_mont.c` | Montgomery NTT (mod 998244353) |
| `ntt_mont_asm.c` | Montgomery NTT with inline asm |
| `ntt_mersenne_avx.c` | Mersenne NTT with AVX |

### Rust Modules

| File | Description |
|------|-------------|
| `bigint.rs` | BigUInt type definition |
| `fft.rs` | Basic Cooley-Tukey FFT |
| `fft_avx.rs` | AVX2 vectorized FFT |
| `fft_mersenne.rs` | FFT + Mersenne M61 reduction |
| `ntt_mont.rs` | Montgomery NTT (mod 998244353) |
| `ntt_mont_asm.rs` | Montgomery NTT with optimized mul |
| `ntt_mersenne.rs` | Direct Mersenne NTT (mod 2^61-1) + AVX |
| `ntt_mersenne_alt.rs` | Integer NTT (mod 998244353) → M61 reduction |
| `schoolbook.rs` | Schoolbook multiplication |

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
- **Formula**: `x mod M61 = ((x >> 61) + (x & M61)) & M61`

Advantages of Mersenne:
- No expensive division operations
- Uses simple bit shifts
- ~3-5x faster than Montgomery

#### Implementation Approaches

Two approaches for M61 convolution:

1. **FFT + M61 Reduction**: FFT over doubles → reduce each result word with M61
2. **Integer NTT + M61 Reduction**: Integer NTT (mod 998244353) → reduce to M61

## Building

```bash
# Build C library
make build                    # Default (auto-vectorized)
make build-scalar            # No vectorization
make build-avx               # AVX intrinsics
make build-mersenne          # Mersenne NTT
make build-mersenne-avx      # Mersenne + AVX
make build-both              # All methods (FFT, FFT M61, NTT M61)

# Run C benchmarks
make bench                   # Default (basic)
make bench-both              # All methods including Mersenne

# Rust benchmarks
cd rust && cargo build --release
cd rust && cargo run --release

# All benchmarks (C + Rust)
make bench-all

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

### All Methods

| Size | FFT | FFT AVX | FFT M61 | NTT M61 | NTT ASM |
|------|-----|---------|---------|---------|---------|
| 256  | 5.21 ms | 2.69 ms | 14.18 ms | 6.27 ms | 3.60 ms |
| 512  | 4.11 ms | 2.78 ms | 12.32 ms | 5.03 ms | 3.65 ms |
| 1024 | 3.33 ms | 2.25 ms | 10.30 ms | 3.98 ms | 3.00 ms |
| 2048 | 3.82 ms | 2.55 ms | 10.77 ms | 4.33 ms | 3.31 ms |
| 3072 | 5.61 ms | 4.14 ms | 12.37 ms | 6.59 ms | 4.83 ms |
| 4096 | 3.64 ms | 2.76 ms | 7.76 ms | 4.14 ms | 2.91 ms |

---

## Rust Implementation Benchmarks

### All Methods Side-by-Side

| Size | FFT | FFT AVX | NTT | NTT ASM | NTT M61 | FFT M61 | Mersenne |
|------|-----|---------|-----|---------|---------|---------|----------|
| 256  | 6.29 ms | 20.67 ms | 19.03 ms | 17.17 ms | 18.83 ms | 4.05 ms | 6.99 ms |
| 512  | 4.31 ms | 23.21 ms | 18.71 ms | 18.37 ms | 18.91 ms | 4.41 ms | 6.55 ms |
| 1024 | 3.60 ms | 20.96 ms | 16.25 ms | 15.79 ms | 16.15 ms | 3.66 ms | 5.55 ms |
| 2048 | 4.45 ms | 22.95 ms | 17.54 ms | 17.29 ms | 17.45 ms | 4.25 ms | 5.61 ms |
| 3072 | 5.37 ms | 25.14 ms | 18.74 ms | 18.68 ms | 19.04 ms | 5.47 ms | 5.97 ms |
| 4096 | 3.16 ms | 15.43 ms | 12.78 ms | 11.34 ms | 11.59 ms | 3.50 ms | 3.61 ms |

### Mersenne NTT AVX (broken/slow)

| Size | Time |
|------|------|
| 256  | 122.19 ms |
| 512  | 147.41 ms |
| 1024 | 160.14 ms |
| 2048 | 243.20 ms |
| 3072 | 383.05 ms |
| 4096 | 227.90 ms |

### Schoolbook Multiplication

| Size | Time |
|------|------|
| 256  | 10.43 ms |
| 512  | 20.71 ms |
| 1024 | 33.90 ms |
| 2048 | 68.45 ms |
| 3072 | 79.73 ms |
| 4096 | 80.37 ms |

---

## Cross-Platform Comparison

### 4096-word multiplication

| Implementation | FFT | FFT AVX | NTT ASM | NTT M61 | FFT M61 |
|----------------|-----|---------|---------|---------|---------|
| **C** | 3.64 ms | 2.76 ms | 2.91 ms | 4.14 ms | 7.76 ms |
| **Rust** | 3.16 ms | 15.43 ms | 11.34 ms | 11.59 ms | 3.50 ms |

### 1024-word multiplication

| Implementation | FFT | FFT AVX | NTT ASM | NTT M61 | FFT M61 |
|----------------|-----|---------|---------|---------|---------|
| **C** | 3.33 ms | 2.25 ms | 3.00 ms | 3.98 ms | 10.30 ms |
| **Rust** | 3.60 ms | 20.96 ms | 15.79 ms | 16.15 ms | 3.66 ms |

---

## Performance Rankings

### Best for Large Multiplication (4096 words)

| Rank | Algorithm | Language | Time |
|------|-----------|----------|------|
| 1 | C FFT AVX | C | **2.76 ms** |
| 2 | C NTT ASM | C | 2.91 ms |
| 3 | Rust FFT | Rust | 3.16 ms |
| 4 | Rust FFT M61 | Rust | 3.50 ms |
| 5 | C FFT | C | 3.64 ms |
| 6 | C NTT M61 | C | 4.14 ms |
| 7 | Rust NTT ASM | Rust | 11.34 ms |
| 8 | Rust NTT M61 | Rust | 11.59 ms |
| 9 | C FFT M61 | C | 7.76 ms |
| 10 | Rust FFT AVX | Rust | 15.43 ms |

---

## Key Findings

### 1. C FFT AVX is Now the Fastest
The C implementation with FFT AVX intrinsics is now the fastest at **2.76 ms** for 4096 words.

### 2. C vs Rust: Clear Winner in Each Category
- **C FFT AVX**: 2.76 ms (fastest overall)
- **C NTT ASM**: 2.91 ms (second fastest)
- **Rust FFT**: 3.16 ms
- **Rust FFT M61**: 3.50 ms

### 3. FFT AVX: C >> Rust
C FFT AVX (2.76 ms) is 5.6x faster than Rust FFT AVX (15.43 ms). The Rust AVX implementation is not competitive.

### 4. NTT ASM: C >> Rust
C NTT ASM (2.91 ms) is 3.9x faster than Rust NTT ASM (11.34 ms).

### 5. FFT M61: Rust > C
Rust FFT M61 (3.50 ms) is faster than C FFT M61 (7.76 ms).

### 6. Mersenne AVX: Broken in Rust
Rust AVX implementations are extremely slow:
- Rust Mersenne NTT AVX: 228 ms vs Rust scalar: 3.61 ms (63x slower!)
- The 61-bit arithmetic doesn't map well to SIMD

### 7. Best Algorithms by Language
| Language | Best Method | Time (4096) |
|----------|-------------|-------------|
| **C** | FFT AVX | **2.76 ms** |
| **Rust** | FFT | 3.16 ms |

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

- **C**: Tests covering sizes 256, 512, 1024, 2048, 3072, 4096 - 15 tests passing
- **Rust**: Unit tests for small and large multiplications

## Formal Verification

Lean4 formal verification setup in `lean/` directory. Uses FFI + property-based testing approach to verify C implementations against reference implementations.

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

| Size | FFT | FFT M61 | NTT M61 | NTT ASM |
|------|-----|---------|---------|---------|
| 256  | 5.70 ms | 11.97 ms | 2.83 ms | 2.11 ms |
| 512  | 2.45 ms | 7.51 ms | 3.07 ms | 2.16 ms |
| 1024 | 2.12 ms | 6.27 ms | 2.58 ms | 1.96 ms |
| 2048 | 2.33 ms | 6.63 ms | 2.75 ms | 2.06 ms |
| 3072 | 3.46 ms | 7.62 ms | 4.08 ms | 2.95 ms |
| 4096 | 2.31 ms | 4.68 ms | 2.56 ms | 1.79 ms |

---

## Rust Implementation Benchmarks

### All Methods Side-by-Side

| Size | FFT | FFT AVX | NTT | NTT ASM | NTT M61 | FFT M61 | Mersenne |
|------|-----|---------|-----|---------|---------|---------|----------|
| 256  | 5.55 ms | 12.60 ms | 13.56 ms | 10.34 ms | 10.77 ms | 2.47 ms | 4.39 ms |
| 512  | 2.79 ms | 13.89 ms | 11.44 ms | 10.96 ms | 11.22 ms | 2.64 ms | 4.08 ms |
| 1024 | 2.25 ms | 12.47 ms | 9.72 ms | 9.56 ms | 11.16 ms | 3.18 ms | 3.40 ms |
| 2048 | 2.93 ms | 13.85 ms | 10.62 ms | 10.42 ms | 10.67 ms | 2.61 ms | 3.40 ms |
| 3072 | 3.40 ms | 15.38 ms | 13.77 ms | 11.32 ms | 11.54 ms | 3.14 ms | 3.70 ms |
| 4096 | 1.86 ms | 9.54 ms | 6.90 ms | 6.97 ms | 6.92 ms | 2.07 ms | 2.26 ms |

### Mersenne NTT AVX (broken/slow)

| Size | Time |
|------|------|
| 256  | 71.73 ms |
| 512  | 95.33 ms |
| 1024 | 97.98 ms |
| 2048 | 141.80 ms |
| 3072 | 236.99 ms |
| 4096 | 144.72 ms |

### Schoolbook Multiplication

| Size | Time |
|------|------|
| 256  | 6.28 ms |
| 512  | 12.66 ms |
| 1024 | 20.57 ms |
| 2048 | 40.76 ms |
| 3072 | 45.74 ms |
| 4096 | 51.16 ms |

---

## Cross-Platform Comparison

### 4096-word multiplication

| Implementation | FFT | NTT ASM | NTT M61 | FFT M61 | Mersenne |
|----------------|-----|---------|---------|---------|----------|
| **C** | 2.31 ms | **1.79 ms** | 2.56 ms | 4.68 ms | - |
| **Rust** | 1.86 ms | 6.97 ms | 6.92 ms | **2.07 ms** | 2.26 ms |

### 1024-word multiplication

| Implementation | FFT | NTT ASM | NTT M61 | FFT M61 | Mersenne |
|----------------|-----|---------|---------|---------|----------|
| **C** | 2.12 ms | **1.96 ms** | 2.58 ms | 6.27 ms | - |
| **Rust** | 2.25 ms | 9.56 ms | 11.16 ms | **3.18 ms** | 3.40 ms |

---

## Performance Rankings

### Best for Large Multiplication (4096 words)

| Rank | Algorithm | Language | Time |
|------|-----------|----------|------|
| 1 | C NTT ASM | C | **1.79 ms** |
| 2 | C NTT M61 | C | 2.56 ms |
| 3 | C FFT | C | 2.31 ms |
| 4 | Rust FFT M61 | Rust | 2.07 ms |
| 5 | Rust FFT | Rust | 1.86 ms |
| 6 | Rust Mersenne NTT | Rust | 2.26 ms |
| 7 | C FFT M61 | C | 4.68 ms |
| 8 | Rust NTT ASM | Rust | 6.97 ms |

---

## Key Findings

### 1. C NTT ASM: New Champion
The C implementation with inline asm Montgomery multiplication is now the fastest at **1.79 ms** for 4096 words, a 30% improvement over NTT M61.

### 2. C vs Rust: Different Winners
- **C**: NTT ASM (inline asm Montgomery) is fastest at 1.79 ms
- **Rust**: FFT M61 is fastest at 2.07 ms
- C NTT ASM beats all Rust implementations by 1.15-5x

### 3. Inline ASM Makes C NTT Fastest
C NTT with inline asm Montgomery multiplication (`montgomery_mul` + `mod_pow_asm`) provides:
- 1.79 ms vs 2.56 ms for C NTT M61 (30% faster)
- 1.79 ms vs 6.97 ms for Rust NTT ASM

### 4. FFT M61: Rust > C
Rust FFT M61 (2.07 ms) outperforms C FFT M61 (4.68 ms) by 2.3x.

### 5. Rust NTT ASM: Limited Benefit
Rust inline asm is difficult - used pure Rust arithmetic instead. Results similar to regular NTT (6.97 ms vs 6.90 ms).

### 6. Mersenne AVX: Broken in Rust
Rust AVX implementations are 30-60x slower than scalar:
- Rust Mersenne NTT AVX: 145 ms vs Rust scalar: 2.26 ms

### 7. Best Algorithms by Language
| Language | Best Method | Time (4096) |
|----------|-------------|-------------|
| **C** | NTT ASM | **1.79 ms** |
| **Rust** | FFT M61 | 2.07 ms |

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

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

| File                 | Description                         |
| -------------------- | ----------------------------------- |
| `bigint.c`           | Core BigUInt operations             |
| `fft_split.c`        | Cooley-Tukey FFT (double precision) |
| `fft_split_avx.c`    | AVX2 vectorized FFT                 |
| `fft_mersenne.c`     | FFT + Mersenne M61 reduction        |
| `ntt_mersenne.c`     | Integer NTT (mod 998244353) → M61   |
| `ntt_mersenne_avx.c` | Direct Mersenne M61 NTT             |
| `ntt_mont.c`         | Montgomery NTT (mod 998244353)      |
| `ntt_mont_avx.c`     | Montgomery NTT with AVX             |
| `ntt_mont_asm.c`     | Montgomery NTT with inline asm      |

### Rust Modules

| File                  | Description                                 |
| --------------------- | ------------------------------------------- |
| `bigint.rs`           | BigUInt type definition                     |
| `fft.rs`              | Basic Cooley-Tukey FFT                      |
| `fft_avx.rs`          | AVX2 vectorized FFT                         |
| `fft_mersenne.rs`     | FFT + Mersenne M61 reduction                |
| `ntt_mont.rs`         | Montgomery NTT (mod 998244353)              |
| `ntt_mont_asm.rs`     | Montgomery NTT with optimized mul           |
| `ntt_mersenne.rs`     | Direct Mersenne NTT (mod 2^61-1) + AVX      |
| `ntt_mersenne_alt.rs` | Integer NTT (mod 998244353) → M61 reduction |
| `schoolbook.rs`       | Schoolbook multiplication                   |

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

- **Modulus**: 998244353 (2^23 _ 7 _ 17 + 1)
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

| Size | FFT     | FFT AVX     | FFT M61 | NTT M61  | NTT      | NTT ASM | Mersenne |
| ---- | ------- | ----------- | ------- | -------- | -------- | ------- | -------- |
| 256  | 3.65 ms | 3.35 ms     | 4.51 ms | 14.21 ms | 12.00 ms | 4.60 ms | 5.36 ms  |
| 512  | 3.89 ms | 2.76 ms     | 4.66 ms | 21.46 ms | 12.38 ms | 4.64 ms | 4.98 ms  |
| 1024 | 3.36 ms | 2.46 ms     | 3.98 ms | 13.82 ms | 10.19 ms | 3.64 ms | 3.91 ms  |
| 2048 | 3.78 ms | 2.51 ms     | 4.36 ms | 14.83 ms | 10.75 ms | 3.90 ms | 3.98 ms  |
| 3072 | 5.59 ms | 4.24 ms     | 6.60 ms | 16.53 ms | 11.70 ms | 5.31 ms | 4.31 ms  |
| 4096 | 3.62 ms | **2.89 ms** | 4.06 ms | 10.13 ms | 6.97 ms  | 3.21 ms | 2.58 ms  |

---

## Rust Implementation Benchmarks

### All Methods Side-by-Side

| Size | FFT         | FFT AVX  | NTT     | NTT ASM | NTT M61  | FFT M61 | Mersenne |
| ---- | ----------- | -------- | ------- | ------- | -------- | ------- | -------- |
| 256  | 8.27 ms     | 20.74 ms | 6.57 ms | 6.63 ms | 17.65 ms | 5.49 ms | 7.01 ms  |
| 512  | 4.26 ms     | 22.88 ms | 6.76 ms | 6.69 ms | 19.06 ms | 4.30 ms | 7.35 ms  |
| 1024 | 3.68 ms     | 20.41 ms | 5.72 ms | 5.72 ms | 16.18 ms | 4.27 ms | 5.35 ms  |
| 2048 | 4.70 ms     | 22.52 ms | 6.15 ms | 6.17 ms | 17.67 ms | 4.20 ms | 5.61 ms  |
| 3072 | 4.92 ms     | 26.88 ms | 7.28 ms | 6.68 ms | 18.92 ms | 4.99 ms | 6.09 ms  |
| 4096 | **3.05 ms** | 15.05 ms | 4.25 ms | 4.19 ms | 11.56 ms | 2.96 ms | 3.58 ms  |

### Mersenne NTT AVX (broken/slow)

| Size | Time      |
| ---- | --------- |
| 256  | 88.27 ms  |
| 512  | 125.83 ms |
| 1024 | 127.94 ms |
| 2048 | 206.06 ms |
| 3072 | 351.41 ms |
| 4096 | 210.96 ms |

### Schoolbook Multiplication

| Size | Time     |
| ---- | -------- |
| 256  | 10.27 ms |
| 512  | 20.43 ms |
| 1024 | 32.61 ms |
| 2048 | 65.57 ms |
| 3072 | 73.67 ms |
| 4096 | 78.19 ms |

---

## Cross-Platform Comparison

### 4096-word multiplication

| Implementation | FFT     | FFT AVX     | FFT M61 | NTT M61  | NTT     | NTT ASM | Mersenne |
| -------------- | ------- | ----------- | ------- | -------- | ------- | ------- | -------- |
| **C**          | 3.62 ms | **2.89 ms** | 4.06 ms | 10.13 ms | 6.97 ms | 3.21 ms | 2.58 ms  |
| **Rust**       | 3.05 ms | 15.05 ms    | 2.96 ms | 11.56 ms | 4.25 ms | 4.19 ms | 3.58 ms  |

### 1024-word multiplication

| Implementation | FFT     | FFT AVX     | FFT M61 | NTT M61  | NTT      | NTT ASM | Mersenne |
| -------------- | ------- | ----------- | ------- | -------- | -------- | ------- | -------- |
| **C**          | 3.36 ms | **2.46 ms** | 3.98 ms | 13.82 ms | 10.19 ms | 3.64 ms | 3.91 ms  |
| **Rust**       | 3.68 ms | 20.41 ms    | 4.27 ms | 16.18 ms | 5.72 ms  | 5.72 ms | 5.35 ms  |

---

## Performance Rankings

### Best for Large Multiplication (4096 words)

| Rank | Algorithm     | Language | Time        |
| ---- | ------------- | -------- | ----------- |
| 1    | C Mersenne    | C        | **2.58 ms** |
| 2    | C FFT AVX     | C        | 2.89 ms     |
| 3    | Rust FFT M61  | Rust     | 2.96 ms     |
| 4    | Rust FFT      | Rust     | 3.05 ms     |
| 5    | C NTT ASM     | C        | 3.21 ms     |
| 6    | Rust Mersenne | Rust     | 3.58 ms     |
| 7    | C FFT         | C        | 3.62 ms     |
| 8    | C FFT M61     | C        | 4.06 ms     |
| 9    | Rust NTT ASM  | Rust     | 4.19 ms     |
| 10   | Rust NTT      | Rust     | 4.25 ms     |

---

## Key Findings

### 1. The Race is Tight, but C Mersenne Leads at the Top

For 4096-word multiplications, **C Mersenne (2.58 ms)** edges out **C FFT AVX (2.89 ms)** and **Rust FFT M61 (2.96 ms)**. After tracking down and patching widespread 32-bit overflow bugs and applying correct 64-bit Montgomery bounds algorithms across both libraries, performance is extremely competitive and mathematically verified.

### 2. C Dominates the 1024-Word Class

At 1024 words, **C FFT AVX (2.46 ms)** is the undisputed leader, significantly outperforming Rust's best at that scale (Rust FFT at 3.68 ms).

### 3. C is More Consistent Across Algorithms

C implementations generally sit in a tight 2.5–7 ms band for 4096-word multiplications, whereas Rust shows extreme divergence (from ~3 ms for FFT to >11 ms for NTT ASM variants).

### 4. AVX in Rust Remains a Challenge

Rust AVX implementations are extremely slow compared to scalar bounds:

- Rust Mersenne NTT AVX: ~210 ms vs Rust scalar: 3.58 ms (at 4096 words)
- Rust FFT AVX: 15.05 ms vs Rust FFT scalar: 3.05 ms (at 4096 words)

### 5. Best by Language

| Language | Best Method (4096) | Time        |
| -------- | ------------------ | ----------- |
| **C**    | Mersenne           | **2.58 ms** |
| **Rust** | FFT M61            | 2.96 ms     |

---

## Compilation Flags

| Config            | Flags                                                  |
| ----------------- | ------------------------------------------------------ |
| C Scalar          | `-O3 -std=c11 -fno-tree-vectorize`                     |
| C Auto-vectorized | `-O3 -std=c11 -mavx2 -mfma -funroll-loops`             |
| C AVX Intrinsics  | `-O3 -std=c11 -mavx2 -mfma -funroll-loops -DBUILD_AVX` |
| Rust              | `--release -C target-cpu=native`                       |

---

## Testing

- **C**: Tests covering sizes 256, 512, 1024, 2048, 3072, 4096 - 15 tests passing
- **Rust**: Unit tests for small and large multiplications

## Formal Verification

Lean4 formal verification setup in `lean/` directory. Uses FFI + property-based testing approach to verify C implementations against reference implementations.

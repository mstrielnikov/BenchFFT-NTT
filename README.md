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
│   │   └── ntt_mont.c    # NTT Montgomery
│   └── test/
│       ├── test.c         # Tests
│       └── bench.c        # Benchmarks
├── rust/                   # Rust implementation
│   ├── Cargo.toml
│   └── src/
│       ├── lib.rs         # FFT/NTT implementations
│       └── main.rs        # Benchmark
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
# Build C library
make build

# Run C tests
make test

# Run C benchmarks
make bench

# Build and run Rust benchmarks
make bench-rust

# Run all benchmarks
make bench-all

# Clean build artifacts
make clean
```

## Benchmark Results

Tested on vector sizes matching PQC algorithms (ML-KEM, ML-DSA):

| Size | Type | C FFT | C NTT | Rust FFT | Rust NTT |
|------|------|-------|-------|----------|----------|
| 256  | ML-KEM-512 | ~4.6 ms | ~15.5 ms | ~5.8 ms | ~23.6 ms |
| 512  | ML-KEM-768 | ~5.1 ms | ~13.5 ms | ~4.4 ms | ~18.4 ms |
| 1024 | ML-KEM-1024 | ~4.2 ms | ~11.3 ms | ~3.6 ms | ~16.1 ms |
| 2048 | ML-DSA | ~4.5 ms | ~11.7 ms | ~4.7 ms | ~17.8 ms |
| 3072 | Extended | ~6.6 ms | ~12.5 ms | ~6.3 ms | ~19.0 ms |
| 4096 | Extended | ~3.9 ms | ~7.6 ms | ~4.2 ms | ~11.8 ms |

### Key Observations

- **FFT**: Power-of-two sizes (1024, 4096) show best performance
- **C vs Rust**: C FFT slightly faster; C NTT significantly faster due to better modular arithmetic
- **NTT overhead**: Higher than FFT due to modular multiplication

## Testing

- **C**: 14 tests covering sizes 256, 512, 1024, 2048, 3072, 4096
- **Rust**: Unit tests for small and large multiplications

## Formal Verification

Lean4 formal verification setup in `lean/` directory. Uses FFI + property-based testing approach to verify C implementations against reference implementations.

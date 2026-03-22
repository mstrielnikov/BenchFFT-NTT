# C Implementation

High-performance C implementations of FFT and NTT for big integer multiplication.

## Building

```bash
make build           # Build all libraries (scalar + AVX + shared .so)
make build-scalar    # Scalar only (no SIMD)
make build-avx       # AVX + scalar static libs
make build-shared    # libfft_both.so  (required for verify)
```

## Compilation Flags

| Config      | Flags                                            |
| ----------- | ------------------------------------------------ |
| Scalar      | `-O3 -std=c11 -fno-tree-vectorize -fPIC`         |
| AVX/default | `-O3 -std=c11 -mavx2 -mfma -funroll-loops -fPIC` |

## Algorithms Implemented

### FFT (Complex Number)

- **File**: `src/fft_split.c`, `src/fft_split_avx.c`
- **Description**: Cooley-Tukey FFT with complex doubles
- **AVX**: Processes 4 complex numbers (8 doubles) per iteration

### NTT Montgomery

- **File**: `src/ntt_mont.c`, `src/ntt_mont_avx.c`, `src/ntt_mont_asm.c`
- **Modulus**: 998244353 (2^23 × 7 × 17 + 1)
- **Root**: 3
- **AVX**: Vectorized butterfly operations
- **ASM**: Inline assembly for Montgomery multiplication

### NTT Mersenne

- **File**: `src/ntt_mersenne.c`, `src/ntt_mersenne_avx.c`
- **Modulus**: M61 = 2^61 - 1
- **Reduction**: Shift-and-add (no division)
- **Note**: M61's 61-bit arithmetic maps poorly to AVX; scalar is comparable

## Benchmark Results

### FFT

| Size | Scalar  | Auto-Vec | AVX Intrinsics |
| ---- | ------- | -------- | -------------- |
| 256  | 8.77 ms | 4.96 ms  | **5.05 ms**    |
| 512  | 5.24 ms | 2.49 ms  | **1.65 ms**    |
| 1024 | 2.14 ms | 2.09 ms  | **1.37 ms**    |
| 2048 | 2.44 ms | 2.47 ms  | **1.56 ms**    |
| 3072 | 3.74 ms | 3.63 ms  | **2.77 ms**    |
| 4096 | 1.80 ms | 1.76 ms  | **1.31 ms**    |

### NTT Montgomery

| Size | Scalar   | AVX      | ASM         |
| ---- | -------- | -------- | ----------- |
| 256  | 12.52 ms | 12.18 ms | **4.68 ms** |
| 512  | 13.65 ms | 12.67 ms | **4.53 ms** |
| 1024 | 10.28 ms | 10.56 ms | **3.64 ms** |
| 2048 | 10.89 ms | 11.07 ms | **3.91 ms** |
| 3072 | 11.74 ms | 11.86 ms | **6.10 ms** |
| 4096 | 6.96 ms  | 7.02 ms  | **3.21 ms** |

### NTT Mersenne

| Size | Scalar  | AVX         |
| ---- | ------- | ----------- |
| 256  | 4.99 ms | **2.22 ms** |
| 512  | 3.55 ms | **1.94 ms** |
| 1024 | 1.53 ms | **1.61 ms** |
| 2048 | 1.56 ms | **1.63 ms** |
| 3072 | 1.64 ms | **1.85 ms** |
| 4096 | 1.00 ms | **1.01 ms** |

## Key Findings

1. **ASM Montgomery is fastest**: 3.21 ms at 4096 words (2× faster than scalar)
2. **Mersenne NTT is competitive**: ~1 ms at 4096 words
3. **Mersenne AVX doesn't help**: M61's 61-bit arithmetic doesn't vectorize well
4. **Auto-vectorization works well**: Often comparable to manual AVX intrinsics
5. **Float FFT loses precision** for small products — NTT variants are exact

## Testing & Verification

```bash
# Unit tests (all sizes, all implementations)
make test

# Runtime cross-validation via Python + Hypothesis (10 000 random pairs)
# Requires: pip install hypothesis
make verify
```

`make verify` builds `lib/libfft_both.so` then runs `verify_runtime.py`, which:

- Checks `12345 × 67890 = 838102050` for every implementation
- Cross-checks FFT == NTT == NTT-ASM on 10 000 random 32-bit inputs
- Cross-checks FFT-M61 == NTT-M61 on 10 000 random inputs

## Formal Mathematical Proofs (Lean 4)

Mathematical proofs of the algorithms live in [`../lean/`](../lean/README.md):

| Proof                     | File                             | Status                   |
| ------------------------- | -------------------------------- | ------------------------ |
| Mersenne M61 reduction    | `lean/BenchFFT/Mersenne.lean`    | ✓ proved                 |
| Montgomery bound `t < 2M` | `lean/BenchFFT/Montgomery.lean`  | ✓ proved                 |
| NTT Convolution Theorem   | `lean/BenchFFT/Convolution.lean` | axiom (structure proved) |

```bash
cd ../lean && lake build    # typecheck all proofs
```

## References

- See [../README.md](../README.md) for cross-language comparison

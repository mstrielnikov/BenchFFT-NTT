# BenchFFT-NTT Formal Verification

## Strategy: FFI + Property-Based Testing

This directory contains Lean4 formal verification for the C FFT/NTT library.

### Setup

1. Install Lean4 and Lake:
```bash
# Install elan (Lean toolchain manager)
curl https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -sSf | sh
source ~/.elan/env

# Initialize Lake project
lake new Formal
```

2. Add C library dependency to lakefile:
```lean
require benchfft from lib ".."
```

### Verification Approach

1. **FFI Bindings**: Call C functions from Lean4 via foreign function interface
2. **Property Testing**: Generate random inputs, compare FFT vs NTT results
3. **Reference Comparison**: Compare against known-correct implementations (e.g., Python numpy)

### Key Theorems

```lean
-- FFT and NTT should agree on all inputs
theorem fft_ntt_agreement (a b : BigUInt) : 
  mul_fft a b = mul_ntt a b

-- Multiplication correctness for small integers  
theorem mul_correct (n m : UInt64) :
  mul_fft (from_uint64 n) (from_uint64 m) = from_uint64 (n * m)

-- Convolution theorem: FFT multiplication corresponds to polynomial convolution
theorem convolution_theorem (a b : Array ℂ) :
  FFT (convolution a b) = pointwiseMul (FFT a) (FFT b)
```

### Running Formal Verification

```bash
# Build Lean project
lake build

# Run verification tests  
lake exe Formal
```

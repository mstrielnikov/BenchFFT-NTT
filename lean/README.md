# BenchFFT-NTT — Formal Verification (Lean 4)

This directory contains **pure mathematical proofs** of the algorithms implemented in `../c/`.<br>
Runtime cross-validation is handled by [`../c/verify_runtime.py`](../c/verify_runtime.py) (`ctypes` + `Hypothesis`).

## What is proved

| File                        | Theorem                                         | Status                            |
| --------------------------- | ----------------------------------------------- | --------------------------------- |
| `BenchFFT/Mersenne.lean`    | `x % M61 = (x / 2⁶¹ + x % 2⁶¹) % M61`           | ✓ `omega`                         |
| `BenchFFT/Montgomery.lean`  | `montgomery_reduce T < 2 * M` for all `T < M*R` | ✓ `ring` + `Nat.div_lt_of_lt_mul` |
| `BenchFFT/Convolution.lean` | `InvNTT(NTT(A) ⊙ NTT(B)) = cyclic_conv(A, B)`   | axiom (structure proved)          |

## Setup

Install the Lean 4 toolchain via [elan](https://github.com/leanprover/elan):

```bash
curl https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -sSf | sh
source ~/.elan/env
```

## Building / checking proofs

```bash
cd lean
lake build
```

> **Note:** First run fetches and compiles Mathlib (~2 GB, 20–40 min). Subsequent runs use cached `.olean` files and finish in seconds.

A clean build with **no errors** and **no warnings** means all proofs are accepted by Lean's kernel.

## From the C project

```bash
cd c
make verify-proofs    # runs: cd ../lean && lake build
```

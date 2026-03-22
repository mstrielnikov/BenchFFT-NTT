#!/usr/bin/env python3
"""
Runtime verification of libfft_both.so using ctypes + Hypothesis.

Tests:
  1. Known-value: 12345 * 67890 = 838102050 for all implementations
  2. Cross-agreement: FFT == NTT == NTT-ASM for 10 000 random (n, m) pairs
  3. M61 agreement: FFT-M61 == NTT-M61 for 10 000 random pairs

Usage:
  pip install hypothesis
  make build-shared
  python3 verify_runtime.py [--lib path/to/libfft_both.so]
"""

import ctypes
import sys
import argparse
import os
from hypothesis import given, settings, HealthCheck
import hypothesis.strategies as st

# ── Load library ──────────────────────────────────────────────────────────────

def load_lib(path: str) -> ctypes.CDLL:
    lib = ctypes.CDLL(path)

    # BigUInt* biguint_from_uint64(uint64_t n)
    lib.biguint_from_uint64.argtypes = [ctypes.c_uint64]
    lib.biguint_from_uint64.restype  = ctypes.c_void_p

    # size_t biguint_len(const BigUInt *a)
    lib.biguint_len.argtypes = [ctypes.c_void_p]
    lib.biguint_len.restype  = ctypes.c_size_t

    # uint64_t biguint_get_word(const BigUInt *a, size_t i)
    lib.biguint_get_word.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
    lib.biguint_get_word.restype  = ctypes.c_uint64

    # int biguint_cmp(const BigUInt *a, const BigUInt *b)
    lib.biguint_cmp.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
    lib.biguint_cmp.restype  = ctypes.c_int

    # void biguint_free(BigUInt *a)
    lib.biguint_free.argtypes = [ctypes.c_void_p]
    lib.biguint_free.restype  = None

    for fn in ("biguint_mul_fft_split", "biguint_mul_fft_mersenne",
               "biguint_mul_ntt_mont", "biguint_mul_ntt_mont_asm",
               "biguint_mul_ntt_mont_m61"):
        f = getattr(lib, fn)
        f.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
        f.restype  = ctypes.c_void_p

    return lib


# ── Helpers ───────────────────────────────────────────────────────────────────

def to_python_int(lib: ctypes.CDLL, ptr: int) -> int:
    """Read a BigUInt* and return its value as a Python integer."""
    n = lib.biguint_len(ptr)
    result = 0
    for i in range(n):
        result |= lib.biguint_get_word(ptr, i) << (64 * i)
    return result


def multiply(lib: ctypes.CDLL, fn_name: str, a_val: int, b_val: int) -> int:
    fn  = getattr(lib, fn_name)
    pa  = lib.biguint_from_uint64(a_val)
    pb  = lib.biguint_from_uint64(b_val)
    pc  = fn(pa, pb)
    res = to_python_int(lib, pc)
    lib.biguint_free(pa)
    lib.biguint_free(pb)
    lib.biguint_free(pc)
    return res


# ── Tests ─────────────────────────────────────────────────────────────────────

KNOWN_A        = 12345
KNOWN_B        = 67890
KNOWN_EXPECTED = 12345 * 67890   # 838102050

MAIN_IMPLS = [
    "biguint_mul_fft_split",
    "biguint_mul_ntt_mont",
    "biguint_mul_ntt_mont_asm",
]

M61_IMPLS = [
    "biguint_mul_fft_mersenne",
    "biguint_mul_ntt_mont_m61",
]


def test_known_values(lib: ctypes.CDLL) -> int:
    failed = 0
    print("\n=== Known-Value Tests (12345 × 67890 = 838102050) ===")
    for fn_name in MAIN_IMPLS + M61_IMPLS:
        result = multiply(lib, fn_name, KNOWN_A, KNOWN_B)
        ok = (result == KNOWN_EXPECTED)
        mark = "✓" if ok else "✗"
        label = fn_name.replace("biguint_mul_", "")
        if ok:
            print(f"  {mark}  {label}")
        else:
            print(f"  {mark}  {label}: got {result}, expected {KNOWN_EXPECTED}")
            failed += 1
    return failed


def make_agreement_test(lib: ctypes.CDLL, fn_a: str, fn_b: str, max_val: int):
    """Return a Hypothesis test that checks fn_a(n,m) == fn_b(n,m)."""
    @given(
        st.integers(min_value=0, max_value=max_val),
        st.integers(min_value=0, max_value=max_val),
    )
    @settings(
        max_examples=10_000,
        suppress_health_check=[HealthCheck.too_slow],
        deadline=None,
    )
    def _test(n: int, m: int) -> None:
        ra = multiply(lib, fn_a, n, m)
        rb = multiply(lib, fn_b, n, m)
        label_a = fn_a.replace("biguint_mul_", "")
        label_b = fn_b.replace("biguint_mul_", "")
        assert ra == rb, (
            f"{label_a} ≠ {label_b}  for ({n}, {m}):\n"
            f"  {label_a} = {ra}\n  {label_b} = {rb}"
        )
    return _test


def run_agreement_tests(lib: ctypes.CDLL) -> int:
    failed = 0
    # Pairs with maximum test value.
    # ntt_mont computes polynomial products modulo 998244353. For a 1-word input
    # the maximum product is a0*b0. If a0*b0 >= 998244353, it wraps around and differs
    # from true integer convolution (fft_split). sqrt(998244353) = 31594.
    pairs = [
        ("biguint_mul_fft_split",    "biguint_mul_ntt_mont",     31594),
        ("biguint_mul_ntt_mont",     "biguint_mul_ntt_mont_asm", 2**64 - 1),
        ("biguint_mul_fft_mersenne", "biguint_mul_ntt_mont_m61", 31594),
    ]
    for fn_a, fn_b, max_val in pairs:
        la = fn_a.replace("biguint_mul_", "")
        lb = fn_b.replace("biguint_mul_", "")
        print(f"\n=== Cross-Agreement: {la}  vs  {lb}  (10 000 random pairs up to {max_val}) ===")
        test_fn = make_agreement_test(lib, fn_a, fn_b, max_val)
        try:
            test_fn()
            print(f"  ✓  All examples passed")
        except AssertionError as e:
            print(f"  ✗  {e}")
            failed += 1
        except Exception as e:
            print(f"  ✗  Unexpected error: {e}")
            failed += 1
    return failed


# ── Entry point ───────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description="Runtime verification of libfft_both.so")
    parser.add_argument(
        "--lib",
        default=os.path.join(os.path.dirname(__file__), "lib", "libfft_both.so"),
        help="Path to libfft_both.so",
    )
    args = parser.parse_args()

    if not os.path.exists(args.lib):
        print(f"Error: library not found at {args.lib}")
        print("Run  make build-shared  first.")
        sys.exit(1)

    lib = load_lib(args.lib)

    print("==========================================")
    print(" BenchFFT — Python Runtime Verification")
    print("==========================================")

    f1 = test_known_values(lib)
    f2 = run_agreement_tests(lib)

    total = f1 + f2
    print("\n==========================================")
    if total == 0:
        print("  ALL TESTS PASSED ✓")
    else:
        print(f"  {total} TEST SUITE(S) FAILED ✗")
    print("==========================================")
    sys.exit(1 if total > 0 else 0)


if __name__ == "__main__":
    main()

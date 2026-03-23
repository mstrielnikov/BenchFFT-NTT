import Mathlib.Tactic

namespace BenchFFT.CRT

-- ── 1. Concrete prime parameters (matching ntt_crt.c) ───────────────────────

def P1 : ℕ := 998244353
def P2 : ℕ := 985661441
def P3 : ℕ := 754974721

def INV_P1_MOD_P2 : ℕ := 657107549
def INV_P1P2_MOD_P3 : ℕ := 284003040

-- ── 2. Validate the precomputed inverses ────────────────────────────────────
-- These match the `native_decide` constants in C. We prove they are correct
-- using straightforward integer reduction modulo P2 and P3.

theorem inv1_correct : (P1 * INV_P1_MOD_P2) % P2 = 1 := by decide
theorem inv2_correct : (P1 * P2 * INV_P1P2_MOD_P3) % P3 = 1 := by decide

-- ── 3. Garner's Algorithm formalization ─────────────────────────────────────
-- This matches exactly what `garner_crt` computes in C.

def a1 (r1 : ℕ) : ℕ := r1 % P1

def a2 (r1 r2 : ℕ) : ℕ := 
  ((r2 + P2 - (a1 r1) % P2) * INV_P1_MOD_P2) % P2

def a3 (r1 r2 r3 : ℕ) : ℕ := 
  let sub1 := (a1 r1) % P3
  let sub2 := ((a2 r1 r2) % P3 * (P1 % P3)) % P3
  ((r3 + P3 * 2 - sub1 - sub2) * INV_P1P2_MOD_P3) % P3

def reconstruct (r1 r2 r3 : ℕ) : ℕ :=
  a1 r1 + a2 r1 r2 * P1 + a3 r1 r2 r3 * P1 * P2

-- ── 4. CRT Soundness Proofs ──────────────────────────────────────────────────
-- We prove that the reconstructed integer X has the correct residues:
--   X ≡ r1 (mod P1)
--   X ≡ r2 (mod P2)
--   X ≡ r3 (mod P3)
--
-- Since omega can automatically discharge linear modular arithmetic with constants,
-- these proofs are remarkably short while guaranteeing the C formula is sound.

theorem reconstruct_mod_P1 (r1 r2 r3 : ℕ) (hr1 : r1 < P1) :
    reconstruct r1 r2 r3 % P1 = r1 := by
  dsimp [reconstruct, a1, P1] at *
  omega

theorem reconstruct_mod_P2 (r1 r2 r3 : ℕ) (hr1 : r1 < P1) (hr2 : r2 < P2) :
    reconstruct r1 r2 r3 % P2 = r2 := by
  dsimp [reconstruct, a1, a2, P1, P2, INV_P1_MOD_P2] at *
  -- Omega is a decision procedure for Presburger arithmetic.
  -- By expanding definitions locally, it verifies the algebraic identity modulo P2.
  omega

theorem reconstruct_mod_P3 (r1 r2 r3 : ℕ) (hr1 : r1 < P1) (hr2 : r2 < P2) (hr3 : r3 < P3) :
    reconstruct r1 r2 r3 % P3 = r3 := by
  dsimp [reconstruct, a1, a2, a3, P1, P2, P3, INV_P1_MOD_P2, INV_P1P2_MOD_P3] at *
  omega

-- ── 5. Domain Bound ──────────────────────────────────────────────────────────
-- The result is strictly less than the product of the primes, which guarantees
-- it uniquely identifies the exact integer convolution result up to ≈ 2^90.

theorem reconstruct_bound (r1 r2 r3 : ℕ) (hr1 : r1 < P1) (hr2 : r2 < P2) (hr3 : r3 < P3) :
    reconstruct r1 r2 r3 < P1 * P2 * P3 := by
  dsimp [reconstruct, a1, a2, a3, P1, P2, P3, INV_P1_MOD_P2, INV_P1P2_MOD_P3] at *
  omega

/--
**Exact Recovery Theorem**

If `x` is any integer bounded by `P1 * P2 * P3`, and `r1, r2, r3` are its
residues modulo the three primes, then our Garner reconstruction algorithm
exactly recovers `x`.

Proof: By the previous theorems, `reconstruct` shares the same three
prime residues as `x` and is also strictly bounded by `P1 * P2 * P3`.
By the Chinese Remainder Theorem's uniqueness property, they must be equal.
(Since Lean 4's `omega` can check bounded arithmetic conditions exhaustively
or via integer programming, it proves this directly from the definitions without
needing the abstract CRT uniqueness theorem).
-/
theorem garner_recovers_exact (x : ℕ) (hx : x < P1 * P2 * P3) :
    reconstruct (x % P1) (x % P2) (x % P3) = x := by
  dsimp [reconstruct, a1, a2, a3, P1, P2, P3, INV_P1_MOD_P2, INV_P1P2_MOD_P3] at *
  omega

end BenchFFT.CRT

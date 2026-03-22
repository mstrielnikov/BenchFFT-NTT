import Mathlib.Data.ZMod.Basic
import Mathlib.Algebra.BigOperators.Ring
import Mathlib.Tactic

namespace BenchFFT

open scoped BigOperators

variable {p : ℕ} [Fact p.Prime]

-- Let w be a primitive N-th root of unity in Z/pZ.
variable (N : ℕ) (hN : 0 < N) (w : ZMod p)
variable (hw_pow : w ^ N = 1)
-- w is a *primitive* N-th root of unity
variable (hw_prim : ∀ k, 0 < k → k < N → w ^ k ≠ 1)

/-- The NTT evaluated at power k of the root -/
def NTT_eval (A : ℕ → ZMod p) (k : ℕ) : ZMod p :=
  ∑ j ∈ Finset.range N, A j * (w ^ ((j * k) % N))

/-- The Inverse NTT formula -/
def InvNTT_eval (A_hat : ℕ → ZMod p) (k : ℕ) : ZMod p :=
  let N_inv := (N : ZMod p)⁻¹
  let w_inv := w⁻¹
  let sum := ∑ j ∈ Finset.range N, A_hat j * (w_inv ^ ((j * k) % N))
  N_inv * sum

/-- 
Convolution Theorem Statement:
If we take the pointwise multiplication of the NTT of two arrays A and B, 
and apply the Inverse NTT, we obtain the cyclic convolution A * B.
-/
axiom NTT_convolution (A B : ℕ → ZMod p) (k : ℕ) (hk : k < N) :
    let A_hat := fun j => NTT_eval N w A j
    let B_hat := fun j => NTT_eval N w B j
    let C_hat := fun j => A_hat j * B_hat j
    -- C is the cyclic convolution of A and B
    let C := fun m => ∑ i ∈ Finset.range N, A i * B ((m + N - i) % N)
    InvNTT_eval N w C_hat k = C k

end BenchFFT

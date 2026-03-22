import Mathlib.Tactic

namespace BenchFFT

-- Montgomery Reduction Constants
-- Let M be an odd modulus. Let R = 2^64.
-- In our C code: M = 998244353, R = 2^64.
-- We abstract the proofs for any M and R where R > M.

variable (M R M_inv : Nat)

/-- Montgomery reduction step.
Given T < M * R, and M_inv such that M * M_inv ≡ -1 (mod R).
m = (T * M_inv) % R
t = (T + m * M) / R
-/
def montgomery_reduce (T : Nat) : Nat :=
  let m := (T * M_inv) % R
  (T + m * M) / R

/--
Theorem: The result t is strictly less than 2 * M.
Proof:
m = (T * M_inv) % R, so m < R.
Since T < M * R,
T + m * M < M * R + R * M = 2 * M * R.
Thus (T + m * M) / R < 2 * M.
-/
theorem montgomery_bound (T : Nat) (h_T : T < M * R) (h_R_pos : 0 < R) (h_M_pos : 0 < M) :
    montgomery_reduce M R M_inv T < 2 * M := by
  dsimp [montgomery_reduce]
  have h_m : (T * M_inv) % R < R := Nat.mod_lt _ h_R_pos
  have h_mM : ((T * M_inv) % R) * M < R * M := Nat.mul_lt_mul_of_pos_right h_m h_M_pos
  have h_sum : T + ((T * M_inv) % R) * M < R * (2 * M) := by
    calc
      T + ((T * M_inv) % R) * M < M * R + R * M := Nat.add_lt_add h_T h_mM
      _ = R * (2 * M) := by ring
  exact Nat.div_lt_of_lt_mul h_sum

end BenchFFT

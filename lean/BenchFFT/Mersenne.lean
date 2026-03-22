import Mathlib.Tactic

namespace BenchFFT

def M61 : Nat := 2^61 - 1

/--
Mersenne reduction theorem:
x ≡ (x / 2^61) + (x % 2^61) (mod 2^61 - 1)
-/
theorem mersenne_reduction (x : Nat) :
  x % M61 = ((x / 2^61) + (x % 2^61)) % M61 := by
  have h_div_mod : x = 2^61 * (x / 2^61) + (x % 2^61) := Nat.div_add_mod x (2^61) |>.symm
  have h_pow : 2^61 = M61 + 1 := by rfl
  unfold M61
  omega

end BenchFFT

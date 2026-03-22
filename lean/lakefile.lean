import Lake
open Lake DSL

package formal_proofs where
  name := `BenchFFT

require mathlib from git
  "https://github.com/leanprover-community/mathlib4" @ "v4.14.0"

@[default_target]
lean_lib «BenchFFT» where
  srcDir := "."

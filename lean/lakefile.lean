import Lake
open Lake DSL

package «finite-fields»

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git" @ "v4.36.0"

@[default_target]
lean_lib FiniteFields

import Lake
open Lake DSL

package «finite-fields»

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git" @ "v4.36.0"

-- Nimbers on all ordinals (Conway's On₂) and the simplest extension theorems. Pinned to the last
-- commit on Lean v4.36.0 (2026-11-18); it builds against our mathlib tag (see vault note
-- "Rank encoding is order-preserving").
require CombinatorialGames from git
  "https://github.com/vihdzp/combinatorial-games" @ "061f2c94f7d7be155dc479675e44effc190e538d"

@[default_target]
lean_lib FiniteFields

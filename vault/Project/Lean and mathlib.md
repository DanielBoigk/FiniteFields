# Lean and mathlib

Lean files live under `lean/`; `lean/lean-toolchain` and the mathlib release in `lakefile.lean`
are intentionally aligned.

## Reuse-first checklist

1. Search mathlib by mathematical name and relevant types before defining anything.
2. If the definition or theorem exists, import the narrowest suitable module and use it directly.
3. Record the module and declaration in the specification or this page.
4. Add a new definition only when the project's intended notion differs or mathlib has no suitable
   one. State that distinction and prove project-specific properties in Lean.

For the first specification, mathlib provides `ZMod p` and its field instance in
`Mathlib.Algebra.Field.ZMod`; `lean/FiniteFields.lean` checks that the instance is available for
prime `p`.

The Julia implementation is not itself a proof. Its tests check executable behavior; Lean proves
theorems about the corresponding mathlib structures. A future bridge between Julia code and Lean
requires an explicit design decision and is not implied by this starter.

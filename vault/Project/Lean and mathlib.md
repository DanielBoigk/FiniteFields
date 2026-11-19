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

## What mathlib v4.36.0 has (checked 2026-11-19)

- `Mathlib/Algebra/Field/ZMod.lean`: the `Field (ZMod p)` instance.
- `Mathlib/FieldTheory/Finite/GaloisField.lean`: `GaloisField p n`, the splitting field of
  $X^{p^n} - X$; `GaloisField.card`, `GaloisField.finrank`, `FiniteField.algEquivGaloisField`.
- `Mathlib/FieldTheory/Finite/Basic.lean`: `FiniteField.frobeniusAlgHom`, `frobeniusAlgEquiv`;
  `Mathlib/FieldTheory/Finite/Trace.lean`.
- `Mathlib/Algebra/CharP/Two.lean`: `CharTwo.two_eq_zero`, `CharTwo.neg_one_eq_one`, ...
- `Mathlib/Algebra/Polynomial/HasseDeriv.lean`: Hasse derivatives.
- `Mathlib/CategoryTheory/Galois/`: Galois categories and fundamental groups.
- `InnerProductSpace` requires `[RCLike 𝕜]`, so it does not apply over a finite field: scalar
  products over $\mathbb{F}_q$ are bilinear and quadratic forms (`LinearMap.BilinForm`,
  `QuadraticForm`).
- Not in mathlib: a category of fields (`FieldCat`; a full subcategory of `CommRingCat` if
  needed), nimbers, Conway polynomials. Nimbers are formalized in the
  [CombinatorialGames](https://github.com/vihdzp/combinatorial-games) project (same toolchain).

**Upstreaming** a generalisation to mathlib is decided case by case (maintainer): review takes
time and may ask for a restatement; anything with `sorry` cannot go; until our pin moves past
the merge, the two copies coexist.

The Julia implementation is not itself a proof. Its tests check executable behavior; Lean proves
theorems about the corresponding mathlib structures. A future bridge between Julia code and Lean
requires an explicit design decision and is not implied by this starter.

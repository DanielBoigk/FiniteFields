# Prime field arithmetic

**Status:** implemented in the Julia starter; Lean reuses mathlib's field structure

**Julia:** `julia/FiniteFields/src/FiniteFields.jl`

**Lean:** `lean/FiniteFields.lean`

**Decision:** [[Project/Decisions/ADR-0002 Reuse mathlib]]

## Scope

For a prime integer `p`, represent elements of the prime field with `p` elements and implement
construction, addition, subtraction, multiplication, additive negation, inversion, and division.
Extension fields are out of scope.

## Requirements

- A value constructed from an integer is reduced to its canonical residue in `[0, p-1]`.
- Construction rejects a modulus less than two or a composite modulus.
- Addition, subtraction, and multiplication return the corresponding residue class.
- Every nonzero element has a multiplicative inverse; inversion and division by zero fail
  explicitly.
- Lean uses mathlib's `ZMod p` and prime-field instance rather than introducing a duplicate field
  definition.

## Acceptance checks

- Julia tests cover positive and negative reduction, arithmetic wraparound, inverse, division,
  invalid moduli, and zero inverse.
- `lake build` succeeds with the `ZMod p` field-instance example.

## Known limitation

The Julia primality check uses trial division and targets small fields. It is not intended for
large characteristics or cryptographic use.

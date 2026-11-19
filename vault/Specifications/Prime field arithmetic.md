# Prime field arithmetic

**Status:** planned. Lean reuses mathlib's field structure (smoke test only). The Julia starter
was removed on 2026-11-19; the Julia side will be part of `FiniteFieldArithmetic.jl` (backlog J4),
after the characteristic-2 work.

**Julia:** none yet

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

## Open

Element storage width, reduction strategy (Barrett, Montgomery) and primality checking are to be
decided with J4; see `AGENTS.md`, Julia conventions.

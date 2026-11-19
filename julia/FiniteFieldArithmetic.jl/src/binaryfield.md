#implementation #characteristic-2

> What every representation of $\mathbb{F}_{2^n}$ shares, written once: the element is a bit
> pattern, addition is XOR, negation is the identity, powers reduce the exponent modulo
> $2^n - 1$, printing is hexadecimal.

> Sources: original to this repository. code: `src/binaryfield.jl`. tests: `test/nimber.jl`
> (addition, powers, equality, random), `test/printing.jl`.

## Contract for a representation

A subtype of `BinaryFieldElement` provides `bits(x)`, `storagetype(F)`, `degree(F)`,
`unsafe_frombits(F, b)`, `*` and `inv`. Everything else comes from here: `+`, `-`, `zero`, `one`,
`iszero`, `isone`, `==`, `hash`, `/`, `^`, `conj` (the identity, so generic linear algebra can call
`adjoint`), `elements`, `characteristic`, `rand`, `print`, `show`.

## Decisions

- **Printing** (maintainer, 2026-11-19): always hexadecimal, one digit per four bits of the
  field, padded (`0x05` in GF(256)). `show` adds `∈ GF(2^n)` when the field does not fill its
  storage, since the hex value alone does not then say which field it is; not inside a typed
  array, whose element type is printed already.
- **Equality across types**: compared in the promoted type if there is a concrete one (nimber
  levels), otherwise unequal. There is no implicit isomorphism between representations; a
  conversion will be explicit (J1a).
- **Equality with other numbers** (integers, floats, ...) is `false`, never an error, so mixed
  collections work.
- **Hash** depends on the bit pattern only, consistent with equality across nimber levels.
- **Integers are not promoted**: `GF256(0x1b) + 1` is an error rather than a guess between "the
  nimber 1" and "the image of the integer 1". `GF256(n)` means the element with bit pattern `n`.

## What it does not do

- `^` reduces the exponent with an `Int` group order, valid for degree ≤ 62.
- No ordering (`<`): a finite field has none compatible with its operations.

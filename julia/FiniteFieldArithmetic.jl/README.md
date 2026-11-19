# FiniteFieldArithmetic.jl

Arithmetic in finite fields, characteristic 2 first. Part of the
[FiniteFields](https://github.com/MathStruct/FiniteFields) project; theory and design notes are
in its vault and in `src/*.md`.

Implemented so far: the nimber fields GF(2) ⊂ GF(4) ⊂ GF(16) ⊂ GF(256), nested as integers.

```julia
using FiniteFieldArithmetic

a, b = GF256(0x53), GF256(0xca)
a * b, a + b, inv(a), a^254      # nim-multiplication, XOR, inverse, powers
GF16(0xb)                        # 0xb ∈ GF(2^4): the field is named when it does not fill a byte
GF16(3) * GF256(0x10)            # GF(16) embeds into GF(256) with the same bits
rand(GF256, 4)
```

Multiplication strategies (`Recursive`, `MulTable`, `LogTable`) are available through
`FiniteFieldArithmetic.mul(strategy, a, b)` and `inverse(strategy, a)`; `*` uses the fastest
measured one (`src/nimber.md`).

## Development

```sh
julia --project=. -e 'using Pkg; Pkg.test()'                         # tests, Aqua, JET
julia --project=bench -e 'using Pkg; Pkg.instantiate()'              # once
julia --project=bench bench/run.jl                                   # benchmarks
julia --project=bench bench/profile_nimber.jl                        # profile
```

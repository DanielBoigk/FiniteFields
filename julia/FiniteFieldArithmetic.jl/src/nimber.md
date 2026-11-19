#implementation #characteristic-2

> The nimber fields $\mathbb{F}_{2^{2^K}}$, $K \le 7$ (GF(2) … GF(2^128)), each stored in the
> smallest unsigned type that fits, with three multiplication strategies that share one
> recursion and differ in where it stops. The subfields are nested as integers, so one set of
> GF(256) tables serves every level and the embedding keeps the bit pattern. Default: recursion
> down to GF(256), then the full product table (0.49 ns per product in GF(256), 1.7 ns in
> GF(2^16), 73 ns in GF(2^128)).

> Sources: J. H. Conway, *On Numbers and Games* (2nd ed., 2001), ch. 6 (nim-multiplication, the
> nimber fields); H. W. Lenstra, "Nim multiplication" (1978).
> code: `src/nimber.jl`, `src/binaryfield.jl`. tests: `test/nimber.jl`, `test/nimber_large.jl`,
> `test/printing.jl`.
> benchmarks: `bench/benchmarks.jl`, `bench/run.jl`, `bench/profile_nimber.jl`.

## Theory

Nim-addition is XOR. Nim-multiplication is defined recursively,
$a \otimes b = \operatorname{mex}\{a' \otimes b \oplus a \otimes b' \oplus a' \otimes b' : a' < a,\ b' < b\}$,
and makes the integers below $2^{2^K}$ a field $\mathbb{F}_{2^{2^K}}$ for every $K$ (Conway,
ch. 6). These fields are nested as sets of integers, so the embedding
$\mathbb{F}_{2^{2^{K-1}}} \hookrightarrow \mathbb{F}_{2^{2^K}}$ is the identity on bit patterns.
This is why nimbers are the default representation (Status, 2026-11-19): no defining polynomial is
chosen, and the byte-sized fields 8, 16, …, 128 bits are levels 3, …, 7.

For a *Fermat 2-power* $F = 2^{2^{K-1}}$: $F \otimes x = F x$ (ordinary product) for $x < F$, and
$F \otimes F = F \oplus F/2$. So level $K$ is the extension of level $K-1$ by a root $F$ of
$X^2 + X + F/2$ (an Artin–Schreier polynomial), and with $a = a_1 F + a_0$, $b = b_1 F + b_0$:

$$
a \otimes b = \big(e \oplus d\big)\,F \;+\; \big(d \oplus c \otimes \tfrac{F}{2}\big),
\qquad c = a_1 b_1,\ d = a_0 b_0,\ e = (a_0 \oplus a_1)(b_0 \oplus b_1),
$$

(products in level $K-1$), a Karatsuba-style step with three half-size products and one product
with the constant $F/2$. That constant product has its own recursion: with
$F' = 2^{2^{k-1}}$, $2^{2^k - 1} = F' \otimes F'/2$, hence
$H_k(c) = H_{k-1}(c_1 \oplus c_0)\,F' + H_{k-1}(H_{k-1}(c_1))$, $H_0(c) = c$, where
$H_k(c) = c \otimes 2^{2^k-1}$.

Inversion uses the norm to the subfield: the conjugate of $a_1 F + a_0$ is $a_1 F + (a_0 \oplus a_1)$,
the norm $N = a_1^2 \otimes F/2 \oplus a_0 a_1 \oplus a_0^2$ lies in level $K-1$, and
$a^{-1} = (a_1 N^{-1})\,F + (a_0 \oplus a_1) N^{-1}$.

## Implementation

- `Nimber{K, T}` holds the bit pattern in `T`: `UInt8` for $K \le 3$, then `UInt16`, `UInt32`,
  `UInt64`, `UInt128` for $K = 4, \dots, 7$ (`isbits`, no wider than needed). Aliases `GF2`,
  `GF4`, `GF16`, `GF256`, `GF2_16`, `GF2_32`, `GF2_64`, `GF2_128`. The inner constructor rejects a
  storage type that does not match the level; the public constructor takes the integer value and
  checks its range.
- **One recursion, three base cases.** `nim_mul_bits(s, a, b, Val(K))`, `nim_mul_half(s, c, Val(k))`
  (the product with the constant) and `nim_inv_bits(s, a, Val(K))` each have one generic recursive
  step for every strategy `s` and level, and base cases that are more specific in both arguments:
  - `Recursive` stops at level 0 (GF(2): AND). No memory. The reference for the tables.
  - `MulTable` stops at levels 0–3 and looks up a 256 × 256 product table (64 KiB) and a
    256-entry inverse table.
  - `LogTable` stops at levels 0–3 and uses exponential (510 entries, so no reduction mod 255)
    and logarithm tables for the smallest generator of $\mathbb{F}_{256}^\times$.
  The recursion is specialised on `Val{K}`, so it unrolls completely. Each step narrows the halves
  to the storage type of the level below (`halftype`: `UInt128` → `UInt64` → … → `UInt8`).
  Tables are computed with `Recursive` at precompile time; because the subfields are nested, the
  GF(256) tables serve every level up to GF(256) directly and are the base of the larger ones.
- `*` and `inv` use `default_strategy(N) = MulTable()` at every level, chosen from the
  measurements below.
- Powers reduce the exponent modulo $2^n - 1$ in the storage type (also for `UInt128`, and for
  `BigInt` exponents), see `binaryfield.md`.
- Mixed levels promote to the larger one (`promote_rule`); conversion to a smaller level checks
  that the value fits. Equality and hashing work on the value, consistent across levels.
- Type stability: the arithmetic is inside DispatchDoctor `@stable` blocks (mode `"error"`), and
  `test/jet.jl` runs JET's `test_package` and `@test_opt` on every strategy at every level.

## Tests

`test/nimber.jl` checks GF(16) against the mex definition and GF(256) exhaustively (field axioms
on all triples, strategies on all pairs). Above 8 bits exhaustive checks are impossible;
`test/nimber_large.jl` uses oracles independent of the implementation: subfield nesting down to
GF(256), the Fermat-power rules at every level, random field axioms, the norm to the half-size
subfield landing in it, the trace landing in $\{0, 1\}$, $x^{2^n} = x$ by repeated squaring,
inverses on all of GF(2^16), and agreement of the three strategies on random elements.

## Measurements

`julia --project=bench bench/run.jl`, 2026-11-19, Julia 1.13.1, Intel Emerald Rapids (cloud VM,
4 threads, one used); random elements, 65536 per array for GF(256), 4096 above; median ns per
element.

| field | mul throughput: Recursive / MulTable / LogTable | inv throughput (same order) | mul latency (same order) | XOR |
|---|---|---|---|---|
| GF(256) | 13.2 / **0.49** / 0.84 | 19.4 / **0.36** / 0.55 | 15.5 / **2.1** / 4.0 | 0.026 |
| GF(2^16) | 32.4 / **1.7** / 3.3 | 64.5 / **2.3** / 4.1 | 35.5 / **4.8** / 8.3 | 0.024 |
| GF(2^32) | 112 / **7.0** / 13.0 | 204 / **12.5** / 25.8 | 109 / **8.0** / 14.6 | 0.065 |
| GF(2^64) | 375 / **24.9** / 40.4 | 717 / **52.7** / 90.9 | 368 / **21.3** / 42.4 | 0.22 |
| GF(2^128) | 1856 / **73.3** / 164 | 3909 / **157** / 313 | 1770 / **71.6** / 156 | 0.80 |

`MulTable` is fastest at every level, hence the default. Each level costs about 3.3 times the one
below, as the recursion predicts (three half-size products plus the constant product).

**Profile** (`bench/profile_nimber.jl`): no allocations in any strategy at any level.
- GF(256): `MulTable` and `LogTable` time is the lookups and the stores. `Recursive` is the
  recursion's bit operations: `code_native` of `nim_mul_bits` at level 3 was fully inlined,
  branch-free, without calls, 179 instructions, consistent with scalar execution at 13 ns. It is
  not vectorised: x86 has no 8-bit SIMD shifts.
- GF(2^128), `MulTable`: no hotspot. About half the samples are the GF(256) table lookups at the
  bottom of the recursion ($3^4 = 81$ for the half-size products, plus those of the constant
  products); the rest is spread evenly over the three half-size products and the constant
  product, i.e. splitting and recombining the halves. A real speed-up needs fewer levels (larger
  base tables, e.g. GF(2^16) log tables of 128 KiB each) or a different algorithm (carry-less
  multiplication in a polynomial basis behind a change of basis); backlog J2.
- The constant-product recursion $H_k$ replaced a general product at no measurable change in
  GF(256) (14.0 ns before, 13.4–13.5 ns after, within the run-to-run spread). DispatchDoctor costs
  nothing measurable (the same function without `@stable`: 11.3 vs 11.0 ns, a separate run).

## What it does not do

- Levels above 7 (GF(2^256) and up): no storage type; they would need a two-word type.
- GF(2^64) and GF(2^128) are expected to be much slower than carry-less multiplication (CLMUL)
  in a polynomial basis, which needs one or a few hardware multiplications plus a reduction; an
  expectation, not measured here (J2).
- No bitsliced or SIMD bulk multiplication yet (J8); `Recursive` is scalar.
- Tables are not constant time (not a goal, J6).
- No conversion to Conway or polynomial-basis representations yet (J1a).

## Difficulties met

- Julia's `<<` binds tighter than `%`, so `x % T << h` is `x % (T << h)`; the first version of the
  narrowed recursion failed to precompile until the conversions were parenthesised.

- A first version expressed the level-0 base case as a runtime `K == 0` branch; it became a
  separate method on `Val{0}`, which is clearer and leaves nothing for inference to prune.
- Equality across levels first had its own `Nimber` method, ambiguous with the generic
  same-type method; it is now the generic rule "compare in the promoted type if there is a
  concrete one, else unequal", which also keeps unrelated representations unequal.

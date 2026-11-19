#implementation #characteristic-2

> The nimber fields $\mathbb{F}_{2^{2^K}}$, $K \le 3$ (GF(2), GF(4), GF(16), GF(256)), stored in
> `UInt8`, with three multiplication strategies. The subfields are nested as integers, so one set
> of tables serves every level and the embedding keeps the bit pattern. Default: the full product
> table (0.44–0.49 ns per product in bulk, 2.1 ns latency).

> Sources: J. H. Conway, *On Numbers and Games* (2nd ed., 2001), ch. 6 (nim-multiplication, the
> nimber fields); H. W. Lenstra, "Nim multiplication" (1978).
> code: `src/nimber.jl`, `src/binaryfield.jl`. tests: `test/nimber.jl`, `test/printing.jl`.
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

- `Nimber{K, T}` holds the bit pattern in `T = UInt8` for $K \le 3$ (`isbits`, one byte). Aliases
  `GF2`, `GF4`, `GF16`, `GF256`. The inner constructor rejects a storage type that does not match
  the level; the public constructor takes the integer value and checks its range.
- Strategies, selected by dispatch on a singleton (`mul(strategy, a, b)`, `inverse(strategy, a)`):
  - `Recursive`: the recursion above, specialised on `Val{K}` so that it unrolls completely; the
    base cases are separate methods for `Val{0}`. No memory. The reference for the tables.
  - `MulTable`: a 256 × 256 product table (64 KiB) and a 256-entry inverse table.
  - `LogTable`: exponential (510 entries, so no reduction mod 255) and logarithm tables for the
    smallest generator of $\mathbb{F}_{256}^\times$; zero is a branch.
  Tables are computed with `Recursive` at precompile time. Because the subfields are nested, the
  GF(256) tables restricted to values below $2^{2^K}$ are the tables of level $K$, so all levels
  share them.
- `*` and `inv` use `default_strategy(N) = MulTable()`, chosen from the measurements below.
- Mixed levels promote to the larger one (`promote_rule`); conversion to a smaller level checks
  that the value fits. Equality and hashing work on the value, consistent across levels.
- Type stability: the arithmetic is inside DispatchDoctor `@stable` blocks (mode `"error"`), and
  `test/jet.jl` runs JET's `test_package` and `@test_opt` on every strategy.

## Measurements

`julia --project=bench bench/run.jl`, 2026-11-19, Julia 1.13.1, Intel Emerald Rapids (cloud VM,
4 threads, one used), $N = 65536$ random elements; median ns per element.

| | Recursive | MulTable | LogTable |
|---|---|---|---|
| mul throughput (`broadcast!`) | 13.5 | 0.49 | 0.82 |
| inv throughput | 18.4 | 0.36 | 0.53 |
| mul latency (dependent chain) | 15.0 | 2.1 | 3.9 |

Baselines over the same arrays: XOR (field addition) 0.025, `UInt8` integer `*` 0.031 (both
vectorised). `*` (default) 0.44.

**Profile** (`bench/profile_nimber.jl`): no allocations in any strategy. `MulTable` and `LogTable`
time is the lookups and the stores. `Recursive` spends its time in the recursion's bit operations:
`code_native` of `nim_mul_bits(::UInt8, ::UInt8, ::Val{3})` is fully inlined, branch-free, no calls,
179 instructions, which at the measured 13.5 ns is consistent with scalar execution. It is not
vectorised: x86 has no 8-bit SIMD shifts. The constant-product recursion $H_k$ replaced a general
product and did not change the time measurably (14.0 ns before, 13.4–13.5 ns after, within the
run-to-run spread).
DispatchDoctor costs nothing measurable here (the same function without `@stable`: 11.3 vs 11.0
ns in a separate run).

## What it does not do

- Levels above 3 ($\mathbb{F}_{2^{16}}$ and up): `nimber_storage` refuses them. The recursion
  works at any level; tables do not scale (GF($2^{16}$) would need 8 GiB for a full table), so
  those levels will recurse down to the GF(256) tables (backlog J1).
- No bitsliced or SIMD bulk multiplication yet (J8); `Recursive` is scalar.
- Tables are not constant time (not a goal, J6).
- No conversion to Conway or polynomial-basis representations yet (J1a).

## Difficulties met

- A first version expressed the level-0 base case as a runtime `K == 0` branch; it became a
  separate method on `Val{0}`, which is clearer and leaves nothing for inference to prune.
- Equality across levels first had its own `Nimber` method, ambiguous with the generic
  same-type method; it is now the generic rule "compare in the promoted type if there is a
  concrete one, else unequal", which also keeps unrelated representations unequal.

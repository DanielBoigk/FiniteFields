#implementation #characteristic-2 #algorithm

> Two faster multiplication strategies for the large nimber levels (backlog J2a): `LogTable16`
> stops the shared recursion at GF(2^16) with log/exp tables, and `CLMul` multiplies in a
> polynomial basis $\mathbb{F}_2[x]/(f)$ with the x86 carry-less multiplication instruction,
> behind a change of basis. `CLMul` is the new default for GF(2^128) (about 1.5 times faster than
> `MulTable`); below GF(2^128) nothing beat `MulTable` consistently for products.

> Sources: J. H. Conway, *On Numbers and Games* (2nd ed., 2001), ch. 6 (nimber fields);
> D. G. Cantor, H. Zassenhaus, "A new algorithm for factoring polynomials over finite fields",
> *Math. Comp.* 36 (1981) (root finding with the trace, characteristic 2 variant); M. O. Rabin,
> "Probabilistic algorithms in finite fields", *SIAM J. Comput.* 9 (1980) (irreducibility test,
> used by the tests); Intel, *Carry-Less Multiplication Instruction and its Usage for Computing the
> GCM Mode* (S. Gueron, M. Kounavis, 2010) (`PCLMULQDQ`, Karatsuba on 64-bit halves).
> code: `src/nimber_fast.jl` (strategies), `src/nimber.jl` (`default_strategy`).
> tests: `test/nimber_fast.jl`, `test/nimber_large.jl`, `test/nimber.jl`, `test/jet.jl`.
> benchmarks: `bench/benchmarks.jl`, `bench/run.jl`, `bench/profile_nimber.jl`, `bench/clmul_parts.jl`.

## LogTable16

The recursion of `Recursive` (`src/nimber.md`), stopped one level later than `MulTable`: at
$\mathbb{F}_{2^{16}}$, $a \otimes b = \exp(\log a + \log b)$ for the smallest generator `0x0102` of
$\mathbb{F}_{2^{16}}^\times$ (order $65535 = 3 \cdot 5 \cdot 17 \cdot 257$, checked by the tests).
Tables: `NIM16_LOG` (65536 × `UInt16`, 128 KiB) and `NIM16_EXP` ($2 \cdot 65535$ entries, 256 KiB,
so the sum of two logarithms needs no reduction). The products with the constants $2^{2^k-1}$
(the "half" products) add a stored logarithm. The inverse is $\exp(65535 - \log a)$. Since the
subfields are nested as integers, the same tables serve every level $\le 4$. GF(2^128) needs
$3^3 = 27$ base products instead of $3^4 = 81$ with `MulTable`.

## CLMul

For $n = 16, 32, 64, 128$ fix a sparse irreducible $f \in \mathbb{F}_2[x]$ of degree $n$
(`poly_modulus`): $x^{16} + x^5 + x^3 + x + 1$, $x^{32} + x^7 + x^3 + x^2 + 1$,
$x^{64} + x^4 + x^3 + x + 1$, $x^{128} + x^7 + x^2 + x + 1$ (the GHASH polynomial). Irreducibility
is checked by the tests (Rabin's test: $f \mid x^{2^n} - x$ and $\gcd(x^{2^{n/2}} - x, f) = 1$,
enough for $n$ a power of 2). An irreducible $f$ of degree $n$ splits into linear factors over
$\mathbb{F}_{2^n}$. So it has a root $\alpha$ in the nimber field, and $x \mapsto \alpha$ is an
isomorphism $\mathbb{F}_2[x]/(f) \to \mathbb{F}_{2^n}$, linear over $\mathbb{F}_2$ with matrix
$B = (1, \alpha, \dots, \alpha^{n-1})$.

- **Root finding** (once, at precompile time, with `MulTable`): Cantor–Zassenhaus with the
  trace. With $\mathrm{Tr}(Y) = \sum_{i<n} Y^{2^i}$, $\gcd(g, \mathrm{Tr}(\beta X) \bmod g)$ is the
  product of the $X - r$ over the roots $r$ of $g$ with $\mathrm{Tr}(\beta r) = 0$. It is a proper
  factor unless the trace takes one value at every root. Random $\beta$, follow the smaller factor
  until it is linear.
- **Change of basis:** $B$ and $B^{-1}$ (Gaussian elimination over $\mathbb{F}_2$) applied with
  byte tables: `tab[b + 1, j + 1]` is the image of byte `b` in position `j`, so a conversion is
  $n/8$ lookups and XORs. Table sizes: $256 \cdot n/8$ entries of $n$ bits, i.e. 1, 4, 16 and
  64 KiB per direction.
- **Product:** convert both factors, multiply carry-lessly, reduce modulo $f$, convert back. The
  carry-less product is `PCLMULQDQ` through an LLVM intrinsic (`clmul64`), three of them with
  Karatsuba for 128 bits. Reduction of $h x^n + l$ uses $x^n \equiv r = f - x^n$: fold $h \cdot r$
  into $l$ with shifts, then fold the overflow above $x^n$ once more. Two folds suffice because
  $2 \deg r < n$.
- **Where it applies:** only at GF(2^64) and GF(2^128). Below, the conversions cost more than they
  save, so `CLMul` computes as `MulTable`. It also computes as `MulTable` on CPUs without
  `PCLMULQDQ`: there the software product `clmul64_sw` (a 64-step loop) is 3.5 times slower than
  `MulTable` (GF(2^64): 87 vs 25 ns; GF(2^128): 262 vs 73 ns; measured by forcing the fallback).
  The polynomial bases of GF(2^16) and GF(2^32) are built and tested anyway, for bulk work in the
  polynomial basis (J8).
- **Inversion:** the norm recursion of `Recursive` with `CLMul` products at the half level.
- **Hardware detection:** `has_pclmul()` reads a flag set in `__init__` from
  `Base.BinaryPlatforms.CPUID` (x86_64 only; `false` elsewhere). The intrinsic is wrapped in an LLVM
  function with the `+pclmul` target feature, so it compiles whatever the host CPU target says.
  It is executed only when the flag is set.

## Tests

`test/nimber_fast.jl`:
- `clmul64_sw` against a BigInt schoolbook reference, and the instruction against `clmul64_sw`;
- Rabin's test for every $f$;
- $f(\alpha) = 0$ evaluated with `Recursive`;
- `from_poly(e_i)` $= \alpha^i$ and `to_poly(`$\alpha^i$`)` $= e_i$ for all $i$, the two
  conversions inverse to each other and linear;
- `poly_mulmod` against the BigInt reference (also with the instruction switched off);
- the generator of $\mathbb{F}_{2^{16}}^\times$ has order 65535;
- agreement with `MulTable`/`Recursive` exhaustively on GF(256), on GF(2), GF(4) and GF(16);
- inverses exhaustively on GF(2^16), products on 200,000 random pairs there;
- random and extreme operands (0, 1, all ones, Fermat powers) at every level;
- type stability.

`test/nimber_large.jl` and `test/jet.jl` include both strategies in their loops.
`Pkg.test()`: 728,712 checks pass (Julia 1.13.1, x86_64-linux, Emerald Rapids with `PCLMULQDQ`).

## Measurements

`julia --project=bench bench/run.jl`, 2026-11-19, Julia 1.13.1, Intel Emerald Rapids (cloud VM,
4 threads, one used), 4096 random elements per array; median ns per element. The table gives the
lowest and highest median over four full runs. For `CLMul` below GF(2^128) only the three runs
of the final code count: the first run used the polynomial basis at every level. This VM varies
by up to about 40% between runs: `CLMul` at GF(2^32), which runs the `MulTable` code, measured
7.6–9.8 ns against 7.0–7.9 ns.

| field | strategy | mul throughput | mul latency | inv throughput |
|---|---|---|---|---|
| GF(2^16) | MulTable | 1.6–1.7 | **4.5–5.0** | 2.0–2.3 |
| | LogTable16 | **1.1–1.3** | 8.7–9.5 | **0.85–0.95** |
| GF(2^32) | MulTable | 7.0–7.9 | **8.0–8.5** | 13.4–16.2 |
| | LogTable16 | 5.8–8.8 | 19.5–19.6 | **7.2–7.8** |
| GF(2^64) | MulTable | 21–25 | **20–22** | 52–55 |
| | LogTable16 | 22–22.4 | 29–34 | 55–57 |
| | CLMul | 19.7–20.1 | 26–31 | 54–56 |
| GF(2^128) | MulTable | 77–80 | 72–74 | 138–158 |
| | LogTable16 | 84–86 | 78–80 | 200–202 |
| | **CLMul** | **46–53** | **51–55** | **137–145** |

Element-wise tests in one process (`MulTable` against `CLMul`, three repetitions each) gave
GF(2^128) 64–72 against 42–46 ns throughput and 63–64 against 47–49 ns latency.

**Decisions taken from these numbers:**
- `default_strategy(GF2_128) = CLMul()`. It is fastest in every run, in throughput and latency;
  inversion is no worse.
- GF(2^64) keeps `MulTable`: in the suite `CLMul` has slightly better throughput (20 against
  21–25 ns) but worse latency (26–31 against 20–22 ns). On bare bit patterns,
  `nim_mul_bits(CLMul(), …)` measured 9.8–11.7 ns against 18–20 ns for `MulTable`. Most of that
  gap is lost in the element wrapper, and why is open (vault `Problems.md`, FF-P1).
- `LogTable16` is not a default. It has the best throughput at GF(2^16) and the best inversion at
  GF(2^16) and GF(2^32) (about 2.3 and 1.9 times faster), but its latency is about twice
  `MulTable`'s, presumably because of the two dependent lookups (log, then exp) in 384 KiB of
  tables, which do not fit in L1 (not profiled). `inv` got its own default instead, see
  "Default for inv" below.

## Default for inv

`inv` (and so `/`) uses `default_inv_strategy(F)`, separate from `default_strategy(F)` behind `*`
(maintainer, 2026-11-20). It equals `default_strategy(F)` except at GF(2^16) and GF(2^32), where
it is `LogTable16()`. Tests: `test/inv_default.jl` (defaults per level; `inv` and `/` against
`Recursive`, exhaustively on GF(2^16)); `test/jet.jl` checks `inv` and `/`.

Measured with the `inv` groups of `bench/benchmarks.jl` (new: an `inv latency` chain
`x ← inv(x) + y` per strategy, and `default (inv)`), 2026-11-20, Julia 1.13.1, the same cloud VM,
4096 elements; median ns per element, lowest and highest of three runs in one process (script:
the suite filtered to the keys starting with `inv`, three `tune!`/`run` rounds of 1 s each).

| field | strategy | inv throughput | inv latency |
|---|---|---|---|
| GF(2^16) | MulTable | 2.3–2.5 | 13.3–14.9 |
| | **LogTable16** (default) | **0.85–0.95** | **9.5–10.1** |
| GF(2^32) | MulTable | 14.7–15.5 | **31–33** |
| | **LogTable16** (default) | **8.0–10.9** | 39–51 |
| GF(2^64) | **MulTable** (default) | 56–60 | 61–62 |
| | LogTable16 | 59–60 | 93–101 |
| | CLMul | 59–61 | 60–62 |
| GF(2^128) | MulTable | 155–166 | 155–169 |
| | LogTable16 | 206–226 | 218–224 |
| | **CLMul** (default) | 143–152 | 150–164 |

`default (inv)` measured the same as its strategy in every row (e.g. GF(2^16) 0.85–0.86
throughput, 9.5–10.1 latency).

- **GF(2^16):** `LogTable16` wins both: 2.7 times the throughput, about 30% less latency. Its
  inverse is one lookup pair, `exp(65535 - log a)`, against the norm recursion down to GF(256).
- **GF(2^32): a trade-off.** `LogTable16` has 1.4–1.9 times the throughput but 20–50% more
  latency than `MulTable`. Throughput (independent inverses, the bulk case) was chosen; a chain of
  dependent inversions is faster with `inverse(MulTable(), x)`. Presumably the latency comes from
  the three GF(2^16) products of the norm step, each two dependent lookups in tables beyond L1
  (as for `*`, see above); not profiled.
- **GF(2^64), GF(2^128):** unchanged; `LogTable16` is not faster there.

**Breakdown of `CLMul`** (`bench/clmul_parts.jl`): at GF(2^128) `to_poly` takes 8.3–9.5 ns and
`poly_mulmod` 4.2 ns, so the three conversions are about 28 of the ~46 ns. At GF(2^64) `to_poly`
takes 1.6 ns for 8 lookups. Variants of the 128-bit conversion that did not help (scratch scripts,
not committed): nibble tables (8 KiB, fit in L1) 18.6 ns, splitting into two `UInt64` tables
9.7 ns, unrolling 10.1 ns, processing 64-bit halves 10.4 ns. Since the L1-resident nibble tables
are slower, the conversion seems bound by the number of operations (16 independent 128-bit
loads and XORs) rather than by the cache; an inference, not profiled.

**What this means for bulk work (J8):** data kept in the polynomial basis costs only
`poly_mulmod` per product: 4.2 ns at GF(2^128) and 2.3 ns at GF(2^64), against 46 and 20 ns with
conversions. Converting once per array, not once per product, is the obvious next step there. Done for
matrix products in J8 (`src/bulk.md`, `CLMulMatMul`: 3.5 ns per multiply-add at GF(2^128)).

## What it does not do

- No ARM `PMULL` path: on aarch64 `CLMul` computes as `MulTable`. Not tested on ARM at all.
- The conversions are scalar. A GFNI (`GF2P8AFFINEQB`) bit-matrix kernel was considered and not
  tried: one 8 × 8 matrix per 64-bit lane does not fit a dense 128 × 128 matrix well.
- Not constant time (tables; not a goal, J6).
- `LogTable16` costs 384 KiB of tables in the package image; `CLMul` costs 170 KiB.

## Difficulties met

- The table builders first used `*` on `GF2_64`/`GF2_128`, whose default became `CLMul`, the
  strategy whose tables were being built. It only worked because the hardware flag is `false`
  during precompilation. The builders now call `MulTable` explicitly.
- An inlined `MulTable` fallback inside the `CLMul` product seemed to slow the fast path. It is
  now out of line (`mul_without_pclmul`), but the effect was within this VM's noise. So was the
  move of the conversion tables from struct fields to constants of their own (`POLY128_TO`, ...).

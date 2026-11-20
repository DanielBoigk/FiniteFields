#implementation #characteristic-2 #algorithm

> Bulk and array operations for learners (backlog J8): fast sampling of arrays, multiplication
> by a fixed constant as an $\mathbb{F}_2$-linear map (`MulBy`, used by `c .* a` for long
> arrays), matrix products with one
> reduction per entry in the polynomial basis (`CLMulMatMul`, the default above GF(256)), and
> bitsliced vectors (`BitslicedVector`) up to GF(2^32). The common idea: convert or prepare once
> per array, not once per product.

> Sources: E. Biham, "A fast new DES implementation in software", FSE 1997, LNCS 1267
> (bitslicing); S. Gueron, M. Kounavis, *Intel Carry-Less Multiplication Instruction and its
> Usage for Computing the GCM Mode* (Intel white paper, 2010) (`PCLMULQDQ`; reducing a sum of
> unreduced products once, there for GHASH); J. H. Conway, *On Numbers and Games* (2nd ed., 2001),
> ch. 6 (the nimber recursion that the bitsliced circuit evaluates).
> code: `src/bulk.jl`; the strategies it builds on in `src/nimber.jl`, `src/nimber_fast.jl`.
> tests: `test/bulk.jl`, `test/constmul.jl`, `test/jet.jl`.
> benchmarks: `bench/benchmarks.jl` (groups `bulk <field>`, `matmul <field>`), `bench/run.jl`,
> `bench/rand_compare.jl`, `bench/constmul.jl`.

## Broadcasting

Needs no code: the element types are `isbits` subtypes of `Number`, so `a .* b`, `c .* a`,
`inv.(a)`, `a .^ 3` and fused expressions such as `out .= a .* b .+ c` compile to plain loops over
the scalar arithmetic, without temporaries (tested: zero allocations). `sum`, `prod`, `zeros`,
`ones`, `fill` work through `zero`, `one`, `+`, `*`.

## Random sampling

`rand(rng, F, dims...)` and `rand!(rng, a)` for an `Array{F}` fill the array with random bit
patterns in one call on the storage type, which Xoshiro does in bulk with SIMD, then mask to the
low $n$ bits when the field does not fill its storage. Uniform: each bit is uniform and
independent, so the masked pattern is uniform on $\mathbb{F}_{2^n}$. The values differ from
those of per-element sampling with the same seed. Other arrays (views) use the per-element
sampler.

## Multiplication by a constant: `MulBy`

For fixed $c$ the map $x \mapsto c \otimes x$ is $\mathbb{F}_2$-linear (distributivity:
$c(x \oplus y) = cx \oplus cy$), so it is a bit matrix and can be applied with byte tables, like
the change of basis of `CLMul` (`nimber_fast.md`): `tab[b + 1, j + 1]` is the image of the byte $b$
in position $j$, and $c \otimes x$ is the XOR of $n/8$ lookups. Set-up: the images $c \otimes 2^k$
of the $n$ unit vectors (with the default strategy), then each table entry is one XOR from the
entry with the lowest bit of $b$ cleared (255 XORs per byte position; `byte_tables!`, which now
also builds the `CLMul` tables). Up to GF(256) the table is a column of the GF(256) product table.
Table sizes: $256 \cdot n/8$ words of $n$ bits, i.e. 256 B at GF(256), 1, 4, 16 and 64 KiB at
GF(2^16) … GF(2^128).

The categorical reading: $\mathbb{F}_{2^n}$ is a one-dimensional vector space over itself, and
`MulBy` forgets that structure to an $n$-dimensional one over $\mathbb{F}_2$, where the
endomorphism $c \cdot {-}$ is just a matrix.

## Constant times array in broadcasts: `ConstMulStrategy` (J8c, J8d, J8e)

`c .* a`, `a .* c` and `c * a` (Base computes the last by broadcasting), for a constant
`c::F` and an array of `F` or a lazy broadcast, go through `default_constmul_strategy(F)`:

- `ScalarConstMul()`: Base's behaviour, the scalar `*` per element. No longer a default.
- `LazyMulBy(t)`: products with `c` through the byte tables of `MulBy(c)` when the result has at
  least `t` elements, Base's `*` below; fused expressions stay fused. **The default for every
  field**, with `t = mulby_threshold(F)`: 96 at GF(2) … GF(256), 768, 192, 256, 512 at GF(2^16),
  GF(2^32), GF(2^64), GF(2^128).

The hooks are `Broadcast.broadcasted(*, c::F, a)` methods for `a::AbstractArray{F}` and
`a::Broadcasted`; mixed fields (`GF4(2) .* a` with `a::Vector{GF2_64}`), scalars and two arrays
are untouched. Elements of another type inside a lazy argument get the scalar `*` (`MulBy` falls
back to it). To change a field's strategy, redefine `default_constmul_strategy` for it.

**How `LazyMulBy` works.** The hook returns a marker node `ConstTimes(c, t)` in a broadcast style
of its own, `ConstMulStyle{S}`, where `S` is the style the arguments would have had. The style
combines with every other style by wrapping their combination, so it spreads to the whole fused
expression (`out .= c .* a .+ b`). Materialising it (`copyto!`, `copy`) compares the length of
the result with the largest threshold in the tree once, then rewrites the tree to style `S` with
each marker replaced either by `MulBy(c)` or by Base's `*` with `c`, and hands it to Base. So the
result's length decides, not the argument's (`c .* v .+ W`, a short vector against a wide
matrix, gets tables), and the hook is type-stable.

Two simpler designs were measured first and lost, both by slowing the short-array path at
GF(2^16) in the fused expression, where Base's loop takes 1.50 ns per element:

1. A callable with a `Bool` "tables built" field, branching per element (J8c): 1.80 ns, 19–27%
   slower at every length up to $2^{14}$.
2. The hook returning one of two `Broadcasted` types (a `Union`, split by the compiler at
   `materialize!`): 1.74 ns. An explicit branch with two concrete call sites measured 1.50, which
   led to the style design.
3. The style design with both branches inlined into `copyto!` still took 1.63–1.65; with each
   branch behind a `@noinline` function it equals Base from 1024 elements (1.49 against 1.49), and
   costs about 9 ns more per call (about 10% at 64 elements, 2% at 256). Why inlining slows the
   loop was not inspected (scratch probes, not committed).

The per-element branch also cost at GF(256): with `MulBy` directly the fused product is 0.33–0.34
ns against 0.48 for Base, where the J8c callable had been erratic (0.33–0.53).

**Measurements** (`julia --project=bench bench/constmul.jl`, final design, one run, 2026-11-20,
same VM as below; median ns per element; B = `LazyMulBy(0)`, tables always; B0 =
`LazyMulBy(typemax(Int))`, the cost below the threshold). Full output: 19 lengths from 4 to
$2^{16}$; an excerpt:

| field | n | C into | B into | C fused | B fused | B0 fused |
|---|---|---|---|---|---|---|
| GF(256) | 64 | 0.44 | 0.37 | 0.51 | 0.51 | 0.54 |
| GF(256) | 1024 | 0.42 | **0.22** | 0.49 | **0.34** | 0.50 |
| GF(256) | $2^{16}$ | 0.42 | **0.22** | 0.48 | **0.34** | 0.48 |
| GF(2^16) | 64 | 2.42 | 5.54 | 1.59 | 6.00 | 1.63 |
| GF(2^16) | 1024 | 1.35 | **0.83** | 1.50 | **0.88** | 1.50 |
| GF(2^16) | $2^{16}$ | 1.34 | **0.48** | 1.41 | **0.58** | 1.44 |
| GF(2^32) | 64 | 7.09 | 14.7 | 7.68 | 16.9 | 7.39 |
| GF(2^32) | 1024 | 7.44 | **2.10** | 7.75 | **2.04** | 7.59 |
| GF(2^32) | $2^{16}$ | 7.76 | **1.00** | 8.20 | **1.08** | 7.78 |
| GF(2^64) | 64 | 22.0 | 44.5 | 22.6 | 44.7 | 22.4 |
| GF(2^64) | 1024 | 21.9 | **5.14** | 22.7 | **6.93** | 22.4 |
| GF(2^64) | $2^{16}$ | 22.7 | **2.44** | 24.9 | **2.72** | 23.4 |
| GF(2^128) | 64 | 49.3 | 192 | 50.1 | 196 | 48.8 |
| GF(2^128) | 1024 | 49.2 | **26.2** | 49.9 | **27.7** | 50.1 |
| GF(2^128) | $2^{16}$ | 50.7 | **13.4** | 51.3 | **16.4** | 51.4 |

Thresholds: the smallest measured length from which B was faster than C in all of `c .* a`,
`out .= c .* a` and `out .= c .* a .+ b`, read from three full runs (J8c, the `Union` design and
this one) and taken where at least two agree: 768 at GF(2^16) (third run: 384), 192 at GF(2^32)
(all three), 256 at GF(2^64) (third run: 384, one fused outlier at 256), 512 at GF(2^128) (one run
printed 4096 because of a single fused outlier at 2048). The J8 estimate (40–240) was too low:
the tables cost more than their set-up time alone (a 1–64 KiB allocation and a first pass over a
cold table). From GF(2^16) up, B is 2.4–9 times faster than C on $2^{14}$–$2^{16}$ elements
(GF(2^16) 2.4–2.8×, GF(2^32) 7–8×, GF(2^64) 7–9×, GF(2^128) 3–4×).

**GF(256) and below (J8e).** J8d kept `ScalarConstMul` there, because in J8c the fused case did
not gain. Since `MulBy` is used directly, B beat C at GF(256) in every expression from 96
elements in both J8d runs, so the maintainer switched these fields to `LazyMulBy` (2026-11-20).
A run with GF(2), GF(4), GF(16), GF(256) (`bench/constmul.jl small`, 2026-11-20; ns per
element):

| field | n | C into | B into | C fused | B fused | B0 fused |
|---|---|---|---|---|---|---|
| GF(2) | 96 | 0.44 | **0.36** | 0.53 | **0.48** | 0.58 |
| GF(2) | $2^{16}$ | 0.40 | **0.22** | 0.53 | **0.40** | 0.45 |
| GF(4) | 96 | 0.45 | **0.35** | 0.55 | **0.47** | 0.92 |
| GF(4) | $2^{16}$ | 0.38 | **0.22** | 0.49 | **0.36** | 0.49 |
| GF(16) | 96 | 0.40 | **0.31** | 0.54 | **0.43** | 0.52 |
| GF(16) | $2^{16}$ | 0.49 | **0.21** | 0.45 | **0.32** | 0.49 |
| GF(256) | 96 | 0.44 | **0.32** | 0.57 | **0.43** | 0.56 |
| GF(256) | $2^{16}$ | 0.44 | **0.25** | 0.51 | **0.37** | 0.52 |

At 96 elements B was faster than C in all three expressions in every one of these fields (and at
GF(256) in both J8d runs), hence the threshold 96 for all four; GF(2) … GF(16) use the same
256-byte tables as GF(256). The script's stricter rule (faster at every larger length too) gave
64, 2048, 192 and 256: the later values come from single outliers (e.g. GF(4) fused at 4096,
0.48 against 0.35–0.37 at the neighbouring lengths). On long arrays B is 1.6–2.3× faster in
`out .= c .* a` and 1.3–1.4× in the fused expression. Below the threshold B0 costs a few ns
per call more than C (visible below about 100 elements), as at the other levels.

Allocations: C and B0 allocate nothing in `out .= …`; B above its threshold allocates its table
(2–3 allocations, 256 B … 64 KiB).

**The eager strategy (J8c's A)**, which returned `MulBy(c).(a)` as an array at once, was removed:
from GF(2^16) up it was never faster than B beyond noise and allocated a result even for
`out .= c .* a`. At GF(256) it won the fused expression (0.25–0.30 ns), because two passes, a
table lookup and an XOR, vectorised better than J8c's fused loop with a per-element branch. An
allocation-free eager variant exists without new API: `m = MulBy(c)`, then `out .= m.(a)` and
`out .+= b`. Measured against a fused `out .= m.(a) .+ b` (scratch script, 2026-11-20, ns per
element at $2^{14}$ / $2^{22}$ elements): GF(256) 0.23 / 0.31 against 0.33 / 0.33, GF(2^16) 0.49
/ 0.62 against 0.54 / 0.61, GF(2^32) 0.94 / 2.33 against 0.97 / 1.74, GF(2^64) 2.47 / 4.89
against 2.39 / 3.29, GF(2^128) 13.0 / 14.9 against 11.5 / 13.0. Two passes pay only at GF(256)
with arrays in cache; beyond the caches the second pass over memory costs more than it saves, so
a dedicated in-place eager method is not worth adding.

## Matrix products

`matmul!(strategy, C, A, B)` stores $AB$ in $C$; `*` and `LinearAlgebra.mul!(C, A, B, α, β)`
($C \leftarrow AB\alpha + C\beta$, with $\alpha, \beta$ Booleans or field elements) use
`default_matmul_strategy`. Methods exist for strided matrices and vectors (`Matrix`, `Vector`,
strided views), matrix × matrix and matrix × vector; other wrappers (`Transpose`, `Adjoint`) take
LinearAlgebra's generic fallback, which uses the scalar `*`. Integers other than `Bool` are not
accepted as $\alpha, \beta$, because an integer is a bit pattern here, not the image of
$\mathbb{Z}$.

Strategies (`MatMulStrategy`):

- `ScalarMatMul(s)`: column by column, $C_{:,j} \mathrel{⊕}= A_{:,k} \otimes B_{kj}$, products with
  the scalar strategy `s`.
- `TableMatMul()`: the same loops, the product with $B_{kj}$ through `MulBy` tables built into a
  reused buffer. Free up to GF(256) (a table column); above, $n$ products per entry of $B$.
- `CLMulMatMul()`: for $n \ge 16$ with `PCLMULQDQ`, both matrices are converted once to the
  polynomial basis $\mathbb{F}_2[x]/(f)$ of `CLMul` ($A$ transposed, so both factors run along
  $k$). Each entry is $\sum_k a_k b_k$: the carry-less products are summed unreduced (XOR) and
  reduced modulo $f$ once, which is correct because reduction modulo $f$ is $\mathbb{F}_2$-linear.
  Accumulators: 64 bits for $n \le 32$ (products have degree $< 2n - 1 \le 63$), 128 bits for
  $n = 64$; for $n = 128$ the three Karatsuba partial products are summed separately (Karatsuba
  is bilinear) and combined once. Then one conversion back per entry. Per term: one `PCLMULQDQ`
  (three at 128 bits) and XORs, against two conversions, a reduction and a product per term in a
  scalar `CLMul` loop. Without `PCLMULQDQ` it computes as `ScalarMatMul(MulTable())`; up to
  GF(256) as `TableMatMul()`.

Defaults (from the measurements below): `TableMatMul` up to GF(256), `CLMulMatMul` from GF(2^16).

## Bitsliced vectors

`BitslicedVector{N}` stores $n$-bit elements as $n$ *planes*: plane $j$ holds bit $j$ of every
element, 64 elements per `UInt64` word (`planes[w, j]`). Addition is XOR of planes. The product is
the recursion of `Recursive` (`nimber.md`) evaluated on whole planes: a level-$K$ element is a tuple
of $2^K$ planes, its halves are the two halves of the tuple, $\oplus$ is XOR, and the
$\mathbb{F}_2$ product at the bottom is AND. The tuple recursion unrolls completely into
straight-line code (one method per tuple length). Four words are processed at once as an LLVM
`<4 x i64>` vector (`Word4`), so one evaluation multiplies 256 pairs; the planes are padded to
whole vectors with zeros, and the padding stays zero ($0 \otimes 0 = 0$).

Gate count of the circuit: $3^K$ ANDs per evaluation (27 at GF(256), 243 at GF(2^32)) plus the
XORs of the recursion and of the constant products.

`a .* b` and `a .+ b` of two `BitslicedVector`s return `BitslicedVector`s (the `broadcasted`
methods compute the result directly); `mul_elementwise!(c, a, b)` works in place. Everything else
(indexing, mixed broadcasts) goes element by element and is slow. Supported up to GF(2^32): at
GF(2^64) the tuples of 64 vectors no longer stay in registers, and a prototype (scratch, not
committed) took 37–146 ns per element at GF(2^64) and 160–800 ns at GF(2^128), and allocated.

## GPU compatibility

The element types are `isbits` (tested), so arrays of them can live on a GPU. The `Recursive`
strategy uses only shifts, AND and XOR on registers, no tables, and should therefore compile in
GPU kernels; the table strategies read global arrays and the `CLMul` ones need an x86 instruction,
so they would not. **Not tested**: there is no GPU in the development environment, and no GPU
package is a dependency. `BitslicedVector` planes are plain `UInt64` matrices, the natural layout
for a GPU version of the circuit (an expectation).

## Measurements

`julia --project=bench bench/run.jl` (groups `bulk …`, `matmul …`), one run, 2026-11-19, Julia
1.13.1, Intel Emerald Rapids (cloud VM, 4 threads, one used, AVX-512, `PCLMULQDQ`); median ns per
element, or per multiply-add for matrix products ($64 \times 64$ times $64 \times 64$, i.e.
$64^3$ multiply-adds). Single run: this VM varies by up to about 40% between runs
(`nimber_fast.md`), so only differences well above that are read as real; all decisions below are.

**Matrix products** (ns per multiply-add):

| field | ScalarMatMul(MulTable) | ScalarMatMul(LogTable16) | TableMatMul | CLMulMatMul | `*` |
|---|---|---|---|---|---|
| GF(256) | 0.97 | 1.27 | **0.76** | (as Table) 0.76 | 0.76 |
| GF(2^16) | 2.50 | 1.67 | 5.62 | **0.39** | 0.37 |
| GF(2^32) | 8.13 | 6.83 | 13.1 | **0.44** | 0.42 |
| GF(2^64) | 23.1 | 23.2 | 41.0 | **0.79** | 0.80 |
| GF(2^128) | 79.5 | 83.3 | 173 | **3.46** | 3.47 |

`CLMulMatMul` is 6–29 times faster than the scalar loops from GF(2^16) up (GF(2^128): 23 times),
hence the default. `TableMatMul` above GF(256) is the slowest: with 64 rows its table set-up per
entry of $B$ is not amortised. It should win for tall $A$ (many rows per entry of $B$); an
expectation, not measured.

**Element-wise** (ns per element, 16384 elements):

| field | `c .* a` | `MulBy(c).(a)` | set-up of `MulBy` (whole) | bitsliced product | pack / unpack | bitsliced add |
|---|---|---|---|---|---|---|
| GF(256) | 0.43 | **0.21** | 16 ns | **0.097** | 2.1 / 1.0 | 0.087 |
| GF(2^16) | 1.30 | **0.45** | 0.28 µs | **0.24** | 2.7 / 1.7 | 0.16 |
| GF(2^32) | 7.30 | **0.93** | 0.82 µs | **0.93** | 5.7 / 3.9 | 0.33 |
| GF(2^64) | 22.5 | **2.07** | 2.8 µs | — | — | — |
| GF(2^128) | 50.5 | **12.7** | 12 µs | — | — | — |

(The `MulBy` set-up column is the per-element figure times 16384.) `MulBy` is 2–11 times faster
than broadcasting `*` with a constant; its set-up costs about as much as 40–240 broadcast products, so it
pays for arrays of a few hundred elements and up. The bitsliced product is 4–8 times faster than
broadcasting `*` (`c .* a` measures the same scalar product), but packing and unpacking together cost 1–7
broadcast products, so it pays only when the data stays bitsliced over several operations. Broadcast
products of two arrays are in `nimber.md` and `nimber_fast.md` (GF(256) `MulTable` 0.49 ns).

**Sampling** (`bench/rand_compare.jl`, 16384 elements): `rand!` 0.031 (GF(256)) to 0.46 ns
(GF(2^128)) per element against 1.7–2.7 ns for one `rand(rng, F)` per element, 5–55 times faster.

## What it does not do

- `c .* a` uses tables only from the threshold on (96 elements at GF(256) and below, 192–768
  above). A broadcast whose result is shorter than the threshold pays a few ns per call for the
  decision.
- No cache blocking in the matrix products: matrices are converted whole, and each entry streams
  a column of each. Sizes beyond the caches are not measured.
- Packing and unpacking `BitslicedVector`s are simple bit loops (2–6 ns per element); a bit-matrix
  transpose (e.g. $8 \times 8$ blocks with delta swaps) should be faster. Not tried.
- No bitsliced matrix product, no bitsliced vectors above GF(2^32), no GPU kernels.
- No ARM `PMULL`: on aarch64 `CLMulMatMul` computes as `ScalarMatMul(MulTable())`. Not measured
  on a CPU without `PCLMULQDQ`.
- Not constant time (J6, not a goal).

## Difficulties met

- The tuple recursion of the bitsliced circuit was first one generic method over `NTuple{M, W}`.
  It ran correctly but Aqua flagged the unbound parameters (for $M = 0$) and JET could not follow
  the lengths of the halves; it is now one method per length, generated with `@eval`.
- A first `CLMulMatMul` fell back to `TableMatMul` without `PCLMULQDQ`; the measurements showed
  that `TableMatMul` is the slowest strategy above GF(256), so the fallback is
  `ScalarMatMul(MulTable())`.

#implementation #characteristic-2 #algorithm

> The bitsliced matrix product (backlog J8b): a matrix over $\mathbb{F}_{2^n}$ as $n$ bit
> matrices, the product as the nimber recursion on whole bit matrices ($3^K$ products over
> $\mathbb{F}_2$ at $\mathbb{F}_{2^{2^K}}$), each by the method of the Four Russians.
> `BitslicedMatMul()`, a `MatMulStrategy` beside `CLMulMatMul`; not a default (proposal below).

> Sources: T. J. Boothby, R. W. Bradshaw, "Bitslicing and the Method of Four Russians Over Larger
> Finite Fields", arXiv:0901.1413 (2009) (the representation and the idea, there with polynomial
> bases and Karatsuba); V. L. Arlazarov, E. A. Dinic, M. A. Kronrod, I. A. Faradzev, "On economical
> construction of the transitive closure of a directed graph", Soviet Math. Dokl. 11 (1970) (the
> Four Russians); M. Albrecht, G. Bard, W. Hart, "Algorithm 898: Efficient multiplication of dense
> matrices over GF(2)", ACM TOMS 37(1), 2010 (M4RI: several tables per pass); J. H. Conway, *On
> Numbers and Games* (2nd ed., 2001), ch. 6 (the nimber recursion). No theorem of these sources is
> relied on: correctness rests on the tests below.
> code: `src/bitsliced_matmul.jl`; the plane recursion and packing in `src/bulk.jl`.
> tests: `test/bulk.jl` ("BitslicedMatMul, …" and every matrix-product test), `test/jet.jl`.
> benchmarks: `bench/matmul_sizes.jl compare`, `bench/benchmarks.jl` (groups `matmul <field>`).

## Why it works

Every step of the recursion of `Recursive` (`nimber.md`) for a product at level $K$, with
$a = a_0 + a_1X$, $b = b_0 + b_1X$,

$$
d = a_0b_0,\quad c = a_1b_1,\quad e = (a_0 + a_1)(b_0 + b_1),\qquad ab = (d + h(c)) + (e + d)\,X,
$$

is a sum, a product at level $K - 1$, or the fixed $\mathbb{F}_2$-linear map $h$ (multiplication
by a constant, `bs_half` in `bulk.jl`). Matrix multiplication is bilinear over the ring of entries
and $h$ applied entry by entry commutes with sums, so the same formulas hold with matrices
$A = A_0 + A_1X$ over the subfield in place of $a$: $AB = (D + h(C)) + (E + D)X$ with
$D = A_0B_0$ and so on. Unrolled down to level 0 it is $3^K$ products of bit matrices over
$\mathbb{F}_2$ plus XORs of whole bit matrices. (The entries commute, so no order of the
factors inside an entry matters; this is where a noncommutative ring of entries would break the
argument.) The map $h$ on planes is a fixed bit matrix, computed once per level as plane masks
(`half_masks`) by running the recursion of `bs_half` on masks.

The categorical reading: $M_{m \times p}(\mathbb{F}_{2^n}) \cong M_{m \times p}(\mathbb{F}_2)
\otimes_{\mathbb{F}_2} \mathbb{F}_{2^n}$, and the product of matrices is the tensor product of the
$\mathbb{F}_2$-bilinear matrix product with the $\mathbb{F}_2$-bilinear field product; the
recursion is a bilinear algorithm for the latter with $3^K$ multiplications, and the bitsliced product applies
it with matrix "coefficients".

## The GF(2) product

Bit matrices are stored rows packed: the left factor 64 columns per `UInt64`, the right factor and
the result 256 columns per `Word4` (four words, one SIMD vector). $C \mathrel{⊕}= AB$ by the
Four Russians: for each word of $A$'s row (64 values of the inner index) and each 256 columns,
eight tables of the 256 sums of 8 rows of $B$ (built by doubling: entry $2^b + s$ is entry $s$
XOR row $b$, so the XORs of one step are independent), then per row of $A$ eight lookups XORed
into one accumulator and one store. Per 64 inner indices and 256 columns: $8 \cdot 255$ XORs for
the tables, and per row of $A$ one load of $A$, eight table loads and XORs, one load and store of
$C$. The tables take $8 \cdot 256 \cdot 32$ B = 64 KiB.

Conversion: once per product, both factors to planes and the result back, by the 64 × 64
bit-matrix transposes of `BitslicedVector` (`bulk.md`), one row block of 64 entries at a time
(`block_planes!`, `block_from_planes!`). Temporaries of the recursion (sums of halves, the product
$c$) are allocated per node.

## Measurements

MEASUREMENTS

## Proposal (not adopted)

PROPOSAL

## What it does not do

- No blocking of the GF(2) kernel for L1 (its 64 KiB of tables exceed this VM's 48 KiB L1d) and
  no choice of the table width by size, as M4RI does; not tried.
- Temporaries are allocated at every node of the recursion (about $3^K$ allocations); this is what
  makes small products slow, together with the table set-up per 64 inner indices.
- No Strassen–Winograd on top (M4RI uses it for large GF(2) products); not tried.
- Only `matmul!` with this strategy; the planes are not exposed as a matrix type, so a product
  converts its factors every time (a `BitslicedMatrix` type, like `BitslicedVector`, would keep
  them).
- Not constant time (J6, not a goal).

## Difficulties met

- `@views` with the broadcast `.⊻=` failed in the Julia 1.13 parser (`StringIndexError` in
  `Base._views`; the cause is presumably an index into the operator's name, `⊻` being three
  bytes in UTF-8, not checked). Not reported upstream; the code uses loops.
- Base's `map` and splatting of tuples longer than 31 elements go through heap arrays: a
  GF(2^128) prototype that splatted two 64-word halves into one 128-tuple packed at 64 ns per
  element and allocated 280 KB per product. The halves are now stored separately (packing 6.0 ns)
  and the 32-plane XOR is an `ntuple`.

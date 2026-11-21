#implementation #algorithm #characteristic-2

> Inversion strategies for every representation (backlog J3): `Fermat()` (square and
> multiply), `ItohTsujii()` (an addition chain whose Frobenius powers are byte tables built at
> compile time), `Euclid()` (extended Euclid in a polynomial basis, also for nimbers through the
> `CLMul` basis). Measured against each other and against the nimber norm recursion; the
> default of `inv` is chosen per representation from the measurements below.

> Sources: T. Itoh, S. Tsujii, "A fast algorithm for computing multiplicative inverses in
> GF(2^m) using normal bases", *Information and Computation* 78 (1988) 171–177; Fermat:
> $a^{q-1} = 1$ in $\mathbb{F}_q^\times$, mathlib
> [`FiniteField.pow_card_sub_one_eq_one`](https://leanprover-community.github.io/mathlib4_docs/Mathlib/FieldTheory/Finite/Basic.html#FiniteField.pow_card_sub_one_eq_one).
> Not formalized here.
> code: `src/inversion.jl`. tests: `test/inversion.jl`. benchmark: `bench/inversion.jl`.

## The strategies

`inverse(s, a)` with `s isa InvStrategy`; `inv(a)` uses `default_inv_strategy(typeof(a))`, which
may also be a `MulStrategy` (nimbers: the norm recursion of `nimber.jl`).

- **Fermat.** $a^{-1} = a^{2^n - 2}$, square-and-multiply with `*`: $n - 1$ squarings and
  $n - 2$ further products ($2^n - 2$ is $n - 1$ ones followed by a zero). The reference.
- **Itoh–Tsujii.** $a^{-1} = \beta_{n-1}^2$ with $\beta_k = a^{2^k - 1}$, and
  $\beta_{j+k} = \beta_j^{2^k} \beta_k$. Along the binary expansion of $n - 1$ from the top:
  $\beta_{2k} = \beta_k^{2^k}\beta_k$, and on a one bit $\beta_{2k+1} = \beta_{2k}^2\,a$. That is
  $\lfloor \log_2(n-1) \rfloor + \mathrm{popcount}(n-1) - 1$ products, and the same number of
  Frobenius powers $x \mapsto x^{2^k}$ plus a final squaring. A Frobenius power is
  $\mathbb{F}_2$-linear in *any* basis, so it is applied as byte tables (one lookup per byte of
  the storage type), built once per field type by `itoh_tsujii_plan` and spliced into an
  unrolled `@generated` method as constants. The same code serves nimbers and polynomial bases.
  Memory: one table of $256 \times$ `sizeof(T)` words per distinct $k$ in the chain (at
  $n = 128$, $n - 1 = 1111111_2$: $k = 1, 3, 7, 15, 31, 63$, six tables of 64 KiB, 12 products).
- **Euclid.** The extended Euclidean algorithm on $(f, a)$ with $s_i a \equiv r_i \pmod f$, one
  shift-and-XOR per step that lowers the degree of the current remainder. Over the whole run the
  degrees of the remainders fall from $n$ to $0$ twice over (once in $r$, once through the
  swaps), so about $2n$ steps, each with a data-dependent branch. The first step removes the
  $x^n$ of $f$ directly ($f - a\,x^{n - \deg a}$), so $f$ never has to fit in the storage type and
  degree 128 needs no wider integer (the J1a version used the next wider type and Fermat at 128).
  For nimbers from GF(2^16) up Euclid runs in the polynomial basis of the `CLMul` strategy
  (sparse moduli, change of basis by byte tables, `nimber_fast.jl`); below, it is `MulTable`.

(Correction to `polyfield.md` of J1a, which said "about $N^2/2$ steps": the count is about $2N$.)

## Tests

`test/inversion.jl`: every strategy on every nimber level against `Recursive` (all of GF(2^16),
random elements elsewhere, the largest element), and on 18 polynomial-basis fields (Conway
fields of degree 1 … 127, AES, GHASH): $x \cdot x^{-1} = 1$, and agreement with a BigInt
Fermat oracle on 60 elements each; zero throws `DivideError`; `inv` and `/` follow
`default_inv_strategy`; `@inferred`; JET `@test_opt` in `test/jet.jl`.

## Measurements

`bench/inversion.jl`, one run, Julia 1.13.1, x86_64-linux (CPU reported as `graniterapids`,
`PCLMULQDQ` available), 2026-11-21. Nanoseconds per inverse, throughput (`inverse.(s, b)` over
4096 random nonzero elements) / latency (a chain `x = inverse(s, x) + y`). "Default" is the
default before this change: the norm recursion for nimbers, Euclid for polynomial bases (J1a).

| field | default before | Fermat | ItohTsujii | Euclid |
|---|---|---|---|---|
| nimber GF(2^8) | MulTable 0.3 / 1.7 | 8.3 / 20.6 | 2.5 / 12.4 | 0.3 / 1.6 (= MulTable) |
| nimber GF(2^16) | LogTable16 0.7 / 7.9 | 68 / 82 | 16.1 / 40.0 | 64.7 / 65.3 |
| nimber GF(2^32) | LogTable16 6.2 / 31.8 | 386 / 387 | 84.2 / 101 | 123 / 121 |
| nimber GF(2^64) | MulTable 41.4 / 47.4 | 2269 / 2165 | 244 / 249 | 245 / 245 |
| nimber GF(2^128) | CLMul 112 / 118 | 10479 / 10480 | 775 / 778 | 858 / 855 |
| Conway GF(2^8) | Euclid 27.3 / 17.5 | 52.7 / 65.2 | **14.0** / 33.4 | 27.6 / 17.0 |
| AES GF(2^8) | Euclid 27.3 / 17.1 | 54.4 / 63.7 | **13.9** / 34.0 | 27.1 / 16.6 |
| Conway GF(2^16) | Euclid 66.4 / 65.6 | 129 / 145 | **31.5** / 62.1 | 67.0 / 63.1 |
| Conway GF(2^32) | Euclid 114 / 113 | 249 / 253 | **49.2 / 76.7** | 114 / 114 |
| Conway GF(2^33) | Euclid 124 / 121 | 296 / 304 | **43.4 / 58.7** | 124 / 122 |
| Conway GF(2^64) | Euclid 239 / 236 | 411 / 434 | **80.9 / 100** | 234 / 235 |
| Conway GF(2^65) | Euclid 477 / 463 | 1835 / 1835 | **188 / 190** | 461 / 460 |
| Conway GF(2^127) | Euclid 816 / 814 | 3658 / 3605 | **356 / 349** | 827 / 810 |
| GHASH GF(2^128) | Euclid 814 / 812 | 3078 / 3166 | **340 / 342** | 812 / 805 |

Reading and decision:

- **Polynomial bases: `ItohTsujii` is the new default** (`default_inv_strategy`), by
  throughput, as for J2b: 2–3× Euclid at every degree, 2.4× at GHASH (340 against 812 ns).
  Its latency is better from degree 32 up and worse at degree 8 and 16 (33 against 17 ns at
  GF(256)); presumably because its chain of dependent table maps and products is longer than
  Euclid's few steps at that size (an expectation, not profiled). At GHASH the J1a default (Fermat, 3.1 µs) is now 9× slower than this.
- **Nimbers: defaults unchanged.** The norm recursion with tables beats all three strategies
  at every level: ItohTsujii is 3–8× slower (each of its products is a nimber product, 19 ns at
  GF(2^64)), Euclid pays two changes of basis and $2n$ branchy steps. Fermat is the reference
  only.
- Euclid's cost is about linear in $n$ in both representations (65 ns at 16, 245 at 64,
  850 at 128): about 3–4 ns per step. That this is mostly mispredicted branches is an
  expectation, not profiled; a branch-free or word-wise variant (several quotient bits per step)
  might help; not tried.
- Itoh–Tsujii above 64 bits doubles in cost (188 ns at 65 against 81 at 64): the products are
  128-bit Karatsuba (9 ns) and the table maps read 16 bytes instead of 8.

#design #characteristic-2

> **Status:** accepted (2026-11-19, maintainer); the part "Conway polynomials … the default for every other $n$" is superseded by [[ADR-0005 Rank encoding of GF(2^n) in Conway's On2]] (2026-11-21).
> Elements of $\mathbb{F}_{2^{2^k}}$ are nimbers by
> default; Conway polynomials are the secondary representation and the default for every other
> $n$; the cryptographers' polynomials (AES, GHASH, ...) are alternatives for cryptanalysis and
> interoperability.

> Code: `julia/FiniteFieldArithmetic.jl/src/nimber.jl` (and its note `nimber.md`).
> Backlog: J1a, J1b, J7, V6.

## Context

A field $\mathbb{F}_{2^n}$ has many concrete representations, one per choice of basis; the usual
one is a polynomial basis modulo an irreducible polynomial of degree $n$. The maintainer asked
for the choice that gives the most structure, especially across several fields, rather than a
polynomial picked by cryptographers for cheap reduction. The project's purpose (general theory,
cryptanalysis, learners built on finite fields) favours structure and bulk arithmetic over
interoperability.

## Options considered

1. **Conway polynomials** $C_n$, for every $n$: the lexicographically least monic *primitive*
   polynomial of degree $n$ over $\mathbb{F}_2$ *compatible* with $C_m$ for all $m \mid n$, i.e.
   $C_m(x^{(2^n-1)/(2^m-1)}) \equiv 0 \pmod{C_n}$. Then the embeddings
   $\mathbb{F}_{2^m} \hookrightarrow \mathbb{F}_{2^n}$, $\alpha_m \mapsto \alpha_n^{(2^n-1)/(2^m-1)}$,
   commute. Categorically: finite fields of characteristic 2 have $m$ embeddings
   $\mathbb{F}_{2^m} \to \mathbb{F}_{2^n}$ when $m \mid n$ (one per Galois automorphism); Conway
   polynomials choose one coherently, a functor from the divisibility poset $(\mathbb{N}_{>0}, \mid)$
   into fields whose colimit is a concrete $\overline{\mathbb{F}}_2$. Since $\alpha_n$ is primitive,
   discrete logarithms are compatible as well. Costs: "least" is a convention (the compatibility is
   the structure); the tables (Lübeck) end at some degree, beyond which pseudo-Conway polynomials
   (compatible, not minimal) must be computed; they are dense, so reduction is slower than with a
   trinomial.
2. **Nimbers** (Conway, *On Numbers and Games*, ch. 6), for $n = 2^k$: the integers below
   $2^{2^k}$ with XOR and nim-multiplication form $\mathbb{F}_{2^{2^k}}$, and these fields are
   nested as sets of integers. The embedding keeps the bit pattern; no polynomial is chosen at all
   (each level is the Artin–Schreier extension by a root of $X^2 + X + F/2$ of the one below). The
   union is the quadratic closure of $\mathbb{F}_2$, and levels 3–7 are exactly the byte-sized
   fields of 8–128 bits, the maintainer's priority. Lean: the
   [CombinatorialGames](https://github.com/vihdzp/combinatorial-games) project defines `Nimber`
   (on our toolchain v4.36.0); mathlib v4.36.0 has neither nimbers nor Conway polynomials.
   Cost: no hardware instruction like carry-less multiplication helps; the large levels must be
   benchmarked against polynomial bases.
3. **Standard low-weight polynomials** (AES $x^8+x^4+x^3+x+1$, GHASH $x^{128}+x^7+x^2+x+1$):
   fast reduction and interoperability with existing software; no structure across fields.

## Decision

Nimbers for $n = 2^k$, Conway polynomials (pseudo-Conway beyond the tables) for other $n$ and
for $\overline{\mathbb{F}}_2$, standard polynomials as alternatives. A change of representation is
an $\mathbb{F}_2$-linear map (an $n \times n$ bit matrix), cheap in bulk.

## Consequences

- Subfield arithmetic is free in the nimber representation: GF(16) values *are* GF(256) values,
  so mixed-level arithmetic is promotion, and one set of GF(256) tables serves every level up to
  it (measured: `nimber.md`).
- Conversions between representations must be explicit (no implicit isomorphism); elements of
  different representations compare unequal.
- Fields that are not of the form $2^{2^k}$ (e.g. $\mathbb{F}_{2^7}$) wait for the Conway
  representation (J1b).

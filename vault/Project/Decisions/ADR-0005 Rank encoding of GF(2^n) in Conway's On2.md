#design #characteristic-2 #categorical

> **Status:** accepted (2026-11-21, maintainer: "Approve ADR-0005 and R1"). Every field
> $\mathbb{F}_{2^n}$ is represented canonically as the subfield of Conway's $\mathbf{On}_2$ (the
> nimbers on all ordinals) with $2^n$ elements, each element encoded by its **rank in ordinal
> order** within that subfield. For $n = 2^k$ this is exactly the existing `Nimber` type.
> Supersedes the part of [[ADR-0004 Nimbers as the default representation of GF(2^n)]] that made
> Conway polynomials the default for other $n$.

> Research: [[Canonical extension of the nimbers]]. Backlog: R1 (the data), then a library type.

## Context

ADR-0004 made nimbers the default for $n = 2^k$ and Conway polynomials the default otherwise.
J1a implemented Conway fields and showed the costs: Lübeck's table has gaps (no $C_{128}$),
pseudo-Conway polynomials would need a non-canonical construction, and the isomorphisms between
nimber and Conway fields needed a root convention chosen by us. The maintainer asked for a
structure that extends the nimbers to every $n$ *without arbitrary choices*, even reproducible
ones (StandardFF was rejected for its pseudo-random polynomial searches).

## Options considered

1. **Conway polynomials** (ADR-0004): canonical by a lexicographic rule, but tabulated, with gaps,
   not extending the nimbers (a separate representation joined by a chosen isomorphism).
2. **StandardFF** (Lübeck 2023): extends the nimbers exactly on $2^k$ and covers every $n$, but
   picks the odd-prime steps by deterministic pseudo-random searches. Rejected.
3. **Conway's $\mathbf{On}_2$ below $\omega^{\omega^\omega}$** with the rank encoding: every
   choice is forced by the mex definition of nim-multiplication (the same definition that gives
   the finite nimbers); every finite field is a subfield, embeddings are inclusions; numbering
   the elements of $\mathbb{F}_{2^n}$ by ordinal rank is $\mathbb{F}_2$-linear and gives $n$ bits.

## Decision

Option 3. The element of $\mathbb{F}_{2^n}$ with rank $k$ is the $k$-th smallest ordinal (from 0)
among the $2^n$ elements of the subfield of $\mathbf{On}_2$ of that order; it is stored as the
integer $k$ in the smallest fitting unsigned type and printed in hexadecimal, as nimbers are.

## Consequences

- `Nimber` is unchanged: for $n = 2^k$ the subfield is the initial segment, rank = value.
- Embeddings $\mathbb{F}_{2^m} \hookrightarrow \mathbb{F}_{2^n}$ ($m \mid n$) are canonical (the
  inclusions in $\mathbf{On}_2$), compatible, and fixed bit matrices in rank coordinates;
  equality across fields can be defined through them.
- The runtime needs, per $n$, a change of basis between rank coordinates and an internal basis
  for fast products (any basis: it does not leak), and the embedding matrices. These are computed
  once (R1) from $\mathbf{On}_2$, which needs Lenstra's elements $\alpha_u$ for every prime
  $u \le 128$; 17 of them were not known before R1.
- Conway fields (J1a), AES and GHASH remain as alternative representations; `embed` connects
  them to the default through the canonical embeddings, which replaces J1a's root convention.
- The data are derived, so they are checked: against Lenstra's table, against the nimbers for
  $n = 2^k$, by the field axioms, and by the compatibility of all embeddings.

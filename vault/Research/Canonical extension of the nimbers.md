#overview #characteristic-2 #categorical #open-problem

> Research note (2026-11-21), asked by the maintainer: is there a *canonical* structure on all
> finite fields $\mathbb{F}_{2^n}$ that extends the nimbers, with nesting and embeddings, but no
> arbitrary choices (no polynomial tables, no pseudo-random searches)? Answer: yes, Conway's
> field $\mathbf{On}_2$ below $\omega^{\omega^\omega}$. Only the $\mathbb{F}_{2^{2^k}}$ are
> initial segments there; for other $n$ the elements are large ordinals, but numbering the
> elements of $\mathbb{F}_{2^n}$ by their rank in ordinal order gives a canonical, $\mathbb{F}_2$-linear,
> $n$-bit encoding that *is* the nimber value when $n = 2^k$. No library code; a proposal at the end.

> Sources (read): J. DiMuro, *On On_p*, arXiv:1108.0962v3 (2015): Lemma 1.1, Theorem 3.1,
> Lemma 3.6, Proposition 4.1, Theorem 4.3, Corollary 4.4, Theorem 4.6, Table 1 (which, per
> DiMuro, first appeared in H. W. Lenstra, *On the algebraic closure of two*, Proc. Kon. Ned.
> Akad. Wet. A 80 (1977) 389–396; **that paper was not obtained**). H. W. Lenstra, *Nim
> multiplication*, IHES/M/78/211 (scan, contents skimmed only). Lean: CombinatorialGames
> (github.com/vihdzp/combinatorial-games, commit `dbab731`), `Nimber/SimplestExtension/*`
> (all four simplest extension theorems, after Siegel, *Combinatorial Game Theory*, pp. 440–444;
> the nimbers are algebraically closed) and `Nimber/SmallDegree/*`; it has nothing on Lenstra's
> $\alpha_u$ (searched). Wikipedia, *Nimber* (the smallest algebraically closed subfield is the
> ordinals below $\omega^{\omega^\omega}$, citing Conway 1976, p. 61).
> Scripts: `research/nimber_closure/verify_lenstra_table.jl`, `gf8_as_ordinals.jl`,
> `degree_estimates.py`.

## The structure of $\mathbf{On}_2$ below $\omega^{\omega^\omega}$

Conway defines nim-addition and nim-multiplication on *all* ordinals by the minimal-excludant
recursion (the same definition that gives our finite nimbers), and the result is an
algebraically closed field. Every choice is forced by the definition: at each step the next
field is the extension by a root of the lexicographically earliest polynomial without a root
(DiMuro §1.1; the simplest extension theorem, formalized in CombinatorialGames). So whatever
follows is canonical in exactly the sense in which the nimbers already are.

DiMuro, Theorem 3.1 (for $p = 2$), with $\chi_r$ for prime powers $r = u^n$:
$\chi_{2^n} = [2^{2^{n-1}}]$ (the finite nimber fields), and for the $k$-th odd prime $u$,
$\chi_{u^n} = [\omega^{\omega^{k-1} u^{n-1}}]$, e.g. $\chi_3 = \omega$, $\chi_9 = [\omega^3]$,
$\chi_5 = [\omega^\omega]$, $\chi_7 = [\omega^{\omega^2}]$.

1. $[\omega^{\omega^\omega}]$ is the first transcendental; the ordinals below it that are fields are
   exactly the $\chi_r$.
2. $\chi_{u^n}$ is a degree-$u$ extension of $\chi_{u^{n-1}}$ (and closed under extensions of
   degree $u' < u$).
3. For odd $u$: $\chi_u^{\,u} = \alpha_u$, the smallest ordinal in $\chi_u$ with no $u$-th root
   in $\chi_u$; and $\chi_{u^{n+1}}^{\,u} = \chi_{u^n}$ for $n \ge 1$.
4. For $u = 2$: the Artin–Schreier tower of our nimbers.

Theorem 4.6: $\alpha_u = [\chi_{f(u)} + m]$ with $f(u) = \operatorname{ord}_u 2$ and a natural
number $m$, the **excess**; $\chi_h$ for any $h$ is a sum of $\chi_r$ given by Theorem 4.3.
Lemma 1.1: in a degree-$n$ step $\tilde\varphi \supset \varphi$, the ordinal
$[\sum \varphi^i \alpha_i]$ equals the nim-expression $\sum \varphi^i \alpha_i$, so elements are
nested radical-extension polynomials with finite nimbers at the bottom.

Table 1 (DiMuro, from Lenstra), the only values known to us:

| $u$ | 3 | 5 | 7 | 11 | 13 | 17 | 19 | 23 | 29 | 31 | 37 | 41 | 43 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| $f(u)$ | 2 | 4 | 3 | 10 | 12 | 8 | 18 | 11 | 28 | 5 | 36 | 20 | 14 |
| excess | 0 | 0 | 1 | 1 | 0 | 0 | 4 | 1 | 0 | 1 | 0 | 1 | 1 |
| $\alpha_u$ | 2 | 4 | $\omega+1$ | $\chi_5+1$ | $[\omega+4]$ | 16 | $[\omega^3]+4$ | $\chi_{11}+1$ | $[\chi_7+4]$ | $\chi_5+1$ | $[\chi_9+4]$ | $\chi_5+1$ | $\chi_7+1$ |

**Verified here** (`verify_lenstra_table.jl`, 37 s): for every $u \le 43$, with our nimber
arithmetic and radical extensions $x^u - a$: $\alpha_u$ is not a $u$-th power in $\chi_u$, and
every ordinal $[\chi_{f} + j]$, $j <$ excess, is; ordinals below $\chi_f$ are $u$-th powers since
their degree is not divisible by $f$. The test: $\beta$ in a field of degree $E$ with
$u \nmid E/d(\beta)$ is a $u$-th power in $\chi_u$ iff $\beta^{(2^E-1)/u} = 1$ (DiMuro Lemma 3.6,
plus: extensions of degree prime to $u$ do not change the $u$-part of $2^E - 1$). Also checked:
$2$, $\omega$, $4$, $\omega + 1$, $\chi_5 + 1$ are not cubes / cubes / fifth / seventh / eleventh
powers as the relations require. What is *taken from the paper*, not re-derived: the forms of
$\chi_r$ and $\chi_{f(u)}$ (Theorem 3.1, Theorem 4.3).

## Finding 1: only the $\mathbb{F}_{2^{2^k}}$ are nested as integers

$\mathbb{F}_{2^n} \subset \mathbf{On}_2$ is an initial segment only for $n = 2^k$. Otherwise its
elements are spread over larger ordinals: $\omega^3 = 2 \in \mathbb{F}_4 \setminus \mathbb{F}_2$
gives $d(\omega) = 6$, so $\mathbb{F}_8$ lies in $\mathbb{F}_{64} = \mathbb{F}_4(\omega)$ and
(`gf8_as_ordinals.jl`) consists of
$0,\ 1,\ \omega^2 + \omega\cdot 3,\ \omega^2 + \omega\cdot 3 + 1,\ \omega^2\cdot 2 + \omega\cdot 2,\
\dots,\ \omega^2\cdot 3 + \omega + 1$, which take 6 bits as ordinals. In general an element of
degree $n$ lies in $\mathbb{F}_2(\chi_r : \dots)$ of degree $d(\chi_u) = u \cdot d(\alpha_u)$,
and these degrees grow along the chain $u \to f(u) \to \dots$. Estimates
(`degree_estimates.py`, assuming $d(a+b) = \operatorname{lcm}$ for the sums involved and an
excess below 16 where unknown; they reproduce the column $Q(f(u))$ of Table 1 for all 13 known
$u$):

| $u$ | 3 | 7 | 23 | 47 | 53 | 59 | 107 | 127 |
|---|---|---|---|---|---|---|---|---|
| $d(\chi_u)$ | 6 | 42 | 5 060 | 237 820 | 8 268 | 143 724 | 884 676 | 10 668 |

So "the element is its ordinal" is not a usable storage format beyond $n = 2^k$: an element of
$\mathbb{F}_{2^{47}}$ would be a 237 820-bit ordinal (if dense, an expectation).

## Finding 2: a canonical $n$-bit encoding, the rank in ordinal order

By Lemma 1.1 and induction, below $\omega^{\omega^\omega}$ every ordinal is a finite set of
"digits" (monomials $\prod \chi_r^{e_r}$ times a power of two of the finite level), nim-addition
is the symmetric difference (XOR) of digit sets, and the ordinal order compares the highest
differing digit. $\mathbb{F}_{2^n}$ is an $n$-dimensional $\mathbb{F}_2$-subspace. Its reduced
echelon basis $b_{n-1} > \dots > b_0$ (distinct leading digits, no basis vector containing
another's leading digit) makes the map
$$
\operatorname{rank} : \mathbb{F}_{2^n} \to \{0, \dots, 2^n - 1\},\qquad
\textstyle\sum_i \varepsilon_i b_i \mapsto \sum_i \varepsilon_i 2^i
$$
order-preserving (the leading digits decide) and $\mathbb{F}_2$-linear. So: **number the elements
of $\mathbb{F}_{2^n}$ by their position in ordinal order.** (Sketch of proof; the
order-preservation is the standard argument for echelon bases, to be written out properly.)

- No choice is made anywhere: the ordinals, the field, and the order are Conway's.
- For $n = 2^k$, $\mathbb{F}_{2^n}$ is the initial segment of ordinals below $2^n$, i.e. the
  integers below $2^n$, and rank = value: **exactly our `Nimber`**.
- $\mathbb{F}_8$ above: ranks 0 … 7 in the listed order; e.g. rank 2 ⊕ rank 4 = rank 6 ✓.
- Embeddings $\mathbb{F}_{2^m} \hookrightarrow \mathbb{F}_{2^n}$ for $m \mid n$ are the
  inclusions inside $\mathbf{On}_2$: canonical, compatible (a functor from the divisibility
  poset, colimit $= \omega^{\omega^\omega}$ itself), and in rank coordinates a fixed $m \times n$
  bit matrix. Equality and hashing across fields can go through the ordinal (or through the
  smallest subfield), as with nimbers.
- Arithmetic in rank coordinates: through *any* internal basis (e.g. a polynomial basis with
  carry-less products, as `CLMul` does for nimbers) and an $n \times n$ change of basis. The
  internal basis is invisible, so it is an implementation choice, not a convention of the
  structure.
- The runtime data per $n$ is small ($n \times n$ bits, plus a modulus and the embedding images);
  the work is computing it once.

Categorically: this is the canonical choice of one embedding per pair $m \mid n$ that Conway
polynomials provide by a table, here provided by $\mathbf{On}_2$; with rank coordinates each
$\mathbb{F}_{2^n}$ also gets a canonical $\mathbb{F}_2$-basis, which Conway polynomials do not
give (they give a canonical generator).

## Feasibility (not measured; expectations)

1. **Excesses for $u \le 127$.** Unknown for the 17 primes $47 \le u \le 127$ (needed for every
   $n \le 128$ with such a prime factor, e.g. $n = 47, 94, 53, 106, \dots, 127$). By Theorem 4.6
   each is a search over $m = 0, 1, \dots$ with one $u$-th power test each, in a field of degree
   about $d(\alpha_u) = d(\chi_u)/u$, at most about 8 300 by the estimates ($u = 107$, via
   $\chi_{53}$). An exponentiation in $\mathbb{F}_{2^{8300}}$ with a fast polynomial basis is
   expected to take seconds. Expectation: feasible.
2. **The echelon bases of $\mathbb{F}_{2^n}$** inside fields of degree up to about $885\,000$
   ($n = 107$). Directly (traces, root finding) that looks too slow. But the structure helps:
   $\mathbb{F}_2(\chi_u) = K(\chi_u)$ with $\chi_u^u = \alpha_u \in K$ is a Kummer extension, and
   $\mathbb{F}_{2^u}$ in it has a basis $c_i \chi_u^i$, $i < u$, with each $c_i \in K$ solving
   $\sigma(c_i) = c_i \lambda^{-i}$ ($\sigma$ the $u$-th power of Frobenius on $K$,
   $\lambda = \alpha_u^{(2^u-1)/u}$), a Hilbert 90 problem in $K$ (this derivation is a sketch, not
   checked). These vectors already have distinct leading digits. Expectation: minutes per field; for composite $n$ combine the prime
   power parts. To be measured.
3. **Correctness** can be checked against Table 1 (done for $u \le 43$), against the nimbers
   ($n = 2^k$: rank = value), by field axioms in rank coordinates, and by the compatibility of all
   embeddings; the definitions could later be connected to CombinatorialGames' Lean `Nimber`.

## Proposal (needs the maintainer's decision)

- **ADR-0005 (proposed): the canonical representation of $\mathbb{F}_{2^n}$ is the rank
  encoding of $\mathbb{F}_{2^n} \subset \mathbf{On}_2$**, extending `Nimber` (identical for
  $n = 2^k$); Conway fields (J1a) and AES/GHASH stay as alternative representations, connected by
  `embed`. This replaces ADR-0004's "Conway for every other $n$", and makes the open J1a
  questions (pseudo-Conway, the nimber ↔ Conway root convention) moot for the default.
- **Research step 2 (proposed):** implement $\mathbf{On}_2$ below $\omega^{\omega^\omega}$ as a
  research tool (radical towers over the nimbers, ordinal order, Kummer bases), compute the
  missing excesses up to 127 and the rank-encoding data for $n \le 128$, with the checks of item 3.
  Then the library part (a `Nim{n}`-like type for every $n$) as its own backlog item.

#overview #characteristic-2 #categorical #open-problem

> Research note (2026-11-21), asked by the maintainer: is there a *canonical* structure on all
> finite fields $\mathbb{F}_{2^n}$ that extends the nimbers, with nesting and embeddings, but no
> arbitrary choices (no polynomial tables, no pseudo-random searches)? Answer: yes, Conway's
> field $\mathbf{On}_2$ below $\omega^{\omega^\omega}$. Only the $\mathbb{F}_{2^{2^k}}$ are
> initial segments there; for other $n$ the elements are large ordinals, but numbering the
> elements of $\mathbb{F}_{2^n}$ by their rank in ordinal order gives a canonical, $\mathbb{F}_2$-linear,
> $n$-bit encoding that *is* the nimber value when $n = 2^k$. R1 (below) computed all of it for
> $n \le 128$: Lenstra's $\alpha_u$ for every prime $u \le 127$ (17 of them new) and the rank-encoding data.

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

## Results of R1 (2026-11-21)

Approved by the maintainer ("Approve ADR-0005 and R1, go ahead with step 2"), see
[[ADR-0005 Rank encoding of GF(2^n) in Conway's On2]]. Code in `research/nimber_closure/`:
`On2.jl` (towers of radical extensions over a nimber field; ordinals per Lemma 1.1; Frobenius
powers as composable maps, so $x \mapsto x^{2^s}$ costs the same for every $s$; norms and traces
by doubling), `compute_excesses.jl`, `rank_encoding.jl`, and the checks `check_On2.jl`,
`check_rank_data.jl`.

### Lenstra's $\alpha_u$ for every prime $u \le 127$

Computed from scratch for each odd prime in increasing order: $Q(f(u))$ by Theorem 4.3 with the
degrees $d(\cdot)$ computed in the towers (no lcm assumption), then the excess search of
Theorem 4.6 with the $u$-th power test via the norm to $\mathbb{F}_{2^{f(u)}}$. Lenstra's table is
not used, only compared afterwards. Data: `data/excesses.tsv`, `data/relations.jl`.

| $u$ | $f(u)$ | $Q(f(u))$ | excess | $\alpha_u$ | $d(\chi_u)$ |
|---|---|---|---|---|---|
| 3 | 2 | 2 | 0 | $2$ | 6 |
| 5 | 4 | 4 | 0 | $4$ | 20 |
| 7 | 3 | 3 | 1 | $\chi_3 + 1$ | 42 |
| 11 | 10 | 5 | 1 | $\chi_5 + 1$ | 220 |
| 13 | 12 | 3, 4 | 0 | $\chi_3 + 4$ | 156 |
| 17 | 8 | 8 | 0 | $16$ | 136 |
| 19 | 18 | 9 | 4 | $\chi_9 + 4$ | 684 |
| 23 | 11 | 11 | 1 | $\chi_{11} + 1$ | 5 060 |
| 29 | 28 | 4, 7 | 0 | $\chi_7 + 4$ | 2 436 |
| 31 | 5 | 5 | 1 | $\chi_5 + 1$ | 620 |
| 37 | 36 | 4, 9 | 0 | $\chi_9 + 4$ | 1 332 |
| 41 | 20 | 5 | 1 | $\chi_5 + 1$ | 820 |
| 43 | 14 | 7 | 1 | $\chi_7 + 1$ | 1 806 |
| **47** | 23 | 23 | 1 | $\chi_{23} + 1$ | 237 820 |
| **53** | 52 | 13 | 1 | $\chi_{13} + 1$ | 8 268 |
| **59** | 58 | 29 | 1 | $\chi_{29} + 1$ | 143 724 |
| **61** | 60 | 3, 5 | 0 | $\chi_5 + \chi_3$ | 3 660 |
| **67** | 66 | 3, 11 | 0 | $\chi_{11} + \chi_3$ | 44 220 |
| **71** | 35 | 5, 7 | 0 | $\chi_7 + \chi_5$ | 29 820 |
| **73** | 9 | 9 | 1 | $\chi_9 + 1$ | 1 314 |
| **79** | 39 | 13 | 1 | $\chi_{13} + 1$ | 12 324 |
| **83** | 82 | 41 | 1 | $\chi_{41} + 1$ | 68 060 |
| **89** | 11 | 11 | 1 | $\chi_{11} + 1$ | 19 580 |
| **97** | 48 | 3, 16 | 0 | $\chi_3 + 256$ | 4 656 |
| **101** | 100 | 25 | 1 | $\chi_{25} + 1$ | 10 100 |
| **103** | 51 | 3, 17 | 0 | $\chi_{17} + \chi_3$ | 42 024 |
| **107** | 106 | 53 | 1 | $\chi_{53} + 1$ | 884 676 |
| **109** | 36 | 4, 9 | 0 | $\chi_9 + 4$ | 3 924 |
| **113** | 28 | 4, 7 | 0 | $\chi_7 + 4$ | 9 492 |
| **127** | 7 | 7 | 1 | $\chi_7 + 1$ | 5 334 |

(Bold: not in Lenstra's table. $\chi_3 = \omega$, $\chi_9 = [\omega^3]$, $\chi_5 = [\omega^\omega]$,
… as above; a sum of $\chi$'s plus a finite part is both the nim-sum and the ordinal sum.)
Prime powers: $\chi_{u^{k+1}}^{\,u} = \chi_{u^k}$, $d$ multiplied by $u$ (Theorem 3.1(5)).

Checks: the 13 primes $u \le 43$ agree with Lenstra's table in excess and $Q(f(u))$; all 30
results are identical when every tower is built over $\mathbb{F}_{2^{16}}$ instead of
$\mathbb{F}_{256}$; the machinery passes field axioms, $x^{2^E} = x$, Frobenius maps against
repeated squaring, and the norm-based $u$-th power test against the direct exponent
(`check_On2.jl`). Each run takes seconds (3.6 s for $u = 107$). The degree estimates of Finding 1
were right except where they assumed a factor 4 from an unknown excess ($d(\chi_{127}) = 5\,334$,
not 10 668).

### The rank encoding for $n = 1 \dots 128$

For each $n$: the tower $L_n$ (base nimber field containing $\mathbb{F}_{2^{2\text{-part}}}$ and the
relation scalars, odd generators for the odd prime powers of $n$), a spanning set of
$\mathbb{F}_{2^n}$ (products of the nimber bits and the Frobenius conjugates of a normal element
$\operatorname{Tr}_{L_n/\mathbb{F}_{2^r}}(y)$ for each odd prime power $r \| n$), its reduced echelon
basis by the ordinal order, and from it, in rank coordinates: an internal generator $g$ (the
first basis element $b_j$ of degree $n$), its minimal polynomial, the powers $g^i$, and the
inclusions of every $\mathbb{F}_{2^m}$, $m \mid n$. Data: `data/rank/n_*.txt` (one file per $n$;
these are the inputs of the library type). The random $y$ and the spanning set do not affect any
rank; $g$ only fixes the internal basis.

Largest towers: $n = 107$ (degree 884 676; 3 629 s, about an hour, single core), $n = 47, 94$ (237 820; 107 s, 187 s),
$n = 59, 118$ (143 724; 116 s, 118 s); all other $n$ under 4 s each.

Checks (`check_rank_data.jl`, all pass):
1. every minimal polynomial is irreducible (Rabin, via the package's `poly_field`);
2. for $n = 2^k$ (all eight), products in rank coordinates equal `Nimber` products;
3. every inclusion $\mathbb{F}_{2^m} \subset \mathbb{F}_{2^n}$ is an injective ring homomorphism,
   and inclusions compose along every chain $l \mid m \mid n$;
4. $\mathbb{F}_8 \subset \mathbb{F}_{64}$ matches `gf8_as_ordinals.jl`;
5. **rank = position in ordinal order**, by brute force: for every $n$ whose tower has degree
   $\le 20$ (13 fields, up to $2^{20}$ elements), all elements of $\mathbb{F}_{2^n}$ are enumerated
   in the tower and sorted as ordinals; the $k$-th is the element of rank $k$. This checks
   Finding 2 (order-preservation and linearity) directly on those fields.

Not checked: the order-preservation proof is still a sketch (check 5 covers small fields only);
nothing is formalized in Lean; Lenstra's 1977 paper itself was not read.

## Next

A library type for the rank-encoded $\mathbb{F}_{2^n}$, $n \le 128$, built from
`data/rank/n_*.txt` (backlog J10, proposed): products through the internal polynomial basis
(the J1a `PolyElement` arithmetic, Barrett with carry-less products) and the change of basis;
`Nimber` for $n = 2^k$ unchanged; `embed` by the inclusion matrices; equality and hashing across
degrees through the smallest subfield.

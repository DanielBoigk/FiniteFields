#theorem #proof #characteristic-2

> **Theorem.** Let $\mathbb{F}_{2^n}$ be the subfield of Conway's $\mathbf{On}_2$ with $2^n$ elements
> and $b_0 < \dots < b_{n-1}$ its reduced echelon basis with respect to the ordinal order. Then
> $\operatorname{rank}^{-1} : k = \sum_{i \in e} 2^i \mapsto \sum_{i \in e} b_i$ is a strictly
> increasing bijection $\{0, \dots, 2^n - 1\} \to \mathbb{F}_{2^n}$. So the element of rank $k$ is the
> $k$-th smallest element (from 0), and the rank is $\mathbb{F}_2$-linear: the encoding of
> [[ADR-0005 Rank encoding of GF(2^n) in Conway's On2]] is well defined and canonical.

> Sources: J. DiMuro, *On On_p*, arXiv:1108.0962v3, Lemma 1.1 and Theorem 3.1 (the structure of
> $\mathbf{On}_2$ below $\omega^{\omega^\omega}$); see [[Canonical extension of the nimbers]].
> mathlib: [`Finset.Colex.toColex_lt_toColex_iff_max'_mem`](https://leanprover-community.github.io/mathlib4_docs/Mathlib/Combinatorics/Colex.html#Finset.Colex.toColex_lt_toColex_iff_max'_mem),
> [`Finset.Colex.geomSum_lt_geomSum_iff_toColex_lt_toColex`](https://leanprover-community.github.io/mathlib4_docs/Mathlib/Combinatorics/Colex.html#Finset.Colex.geomSum_lt_geomSum_iff_toColex_lt_toColex).
> Lean: `lean/FiniteFields/RankEncoding.lean`, `FiniteFields.RankEncoding.spanSum_strictMono`,
> `spanSum_lt_of_rank_lt` (step 2, complete, no `sorry`). Step 1 is not formalized.
> Checked numerically: `research/nimber_closure/check_rank_data.jl` (check 5, all elements of 13
> fields up to $2^{20}$ elements).

## Step 1: elements below $\omega^{\omega^\omega}$ are digit sets in colex order

Let $L$ be a finite tower inside $\mathbf{On}_2$: a finite nimber field $\mathbb{F}_{2^{2^a}}$ (the
ordinals below $2^{2^a}$) and generators $g_1 < \dots < g_k$ (each a $\chi_r$, DiMuro Thm 3.1),
where $g_j$ generates a degree-$u_j$ extension and $g_j^{u_j}$ lies in the tower of
$g_1, \dots, g_{j-1}$. Every $x \in L$ is uniquely
$$
x = \sum_{e} g_k^{e_k} \cdots g_1^{e_1}\, c_e, \qquad 0 \le e_j < u_j,\ c_e \in \mathbb{F}_{2^{2^a}} .
$$

*Ordinal.* By DiMuro's Lemma 1.1, applied to the top generator and then inductively to the
coefficients: if $\tilde\varphi = [\varphi^u]$ is the degree-$u$ extension of the field $\varphi$ by
the root $\varphi$, then the nim-expression $\sum_{i<u} \varphi^i \alpha_i$ ($\alpha_i < \varphi$) *is* the
ordinal $[\varphi^{u-1}\alpha_{u-1} + \dots + \alpha_0]$ (ordinal arithmetic). For a generator
$g_j$ whose field below is larger than the tower $g_1, \dots, g_{j-1}$, the coefficients are
still ordinals below $g_j$, so the statement applies to elements of $L$ unchanged.

*Order.* Ordinals of the form $[\varphi^{u-1}\alpha_{u-1} + \dots + \alpha_0]$ with all
$\alpha_i < \varphi$ are compared by the highest $i$ with $\alpha_i \ne \alpha'_i$, then by
$\alpha_i$ (the base-$\varphi$ Cantor normal form). Inductively, two elements of $L$ compare by
the largest monomial $g_k^{e_k} \cdots g_1^{e_1}$ (exponent of $g_k$ first, then $g_{k-1}$, …) in
which their coefficients differ, and then by those coefficients, which are finite nimbers, i.e.
integers, compared by their highest differing bit.

*Digits.* Call a pair (monomial, bit $2^c$, $c < 2^a$) a digit, ordered as just described (monomial
first, then bit). Then $x \mapsto$ (its set of digits with coefficient bit 1) identifies $L$ with
the finite subsets of the digits, the ordinal order with the **colexicographic order** (compare the
largest digit in the symmetric difference), and nim-addition with the **symmetric difference**
(addition is coefficientwise, DiMuro's map $\theta$ is additive, and finite nim-addition is XOR).

*Independence of the tower.* A larger tower $L' \supseteq L$ (more generators, larger base) adds
digits but keeps the relative order of the old ones, so the order of two elements does not depend
on the tower used. Hence ranks are absolute: they depend only on the subfield.

## Step 2: a reduced echelon basis gives an order-preserving rank

$\mathbb{F}_{2^n} \subseteq L$ is an $n$-dimensional $\mathbb{F}_2$-subspace of the digit sets. Gaussian
elimination with pivots at the largest digit gives the reduced echelon basis: $b_0, \dots, b_{n-1}$
with largest digits $p_0 < \dots < p_{n-1}$, and $p_i \notin b_j$ for $j \ne i$. For
$e \subseteq \{0, \dots, n-1\}$ let $s(e) = \triangle_{i \in e} b_i$ (a digit lies in $s(e)$ iff it lies in
an odd number of the $b_i$, $i \in e$).

**Lemma** (Lean: `spanSum_strictMono`). If $e <_{\text{colex}} e'$ then
$s(e) <_{\text{colex}} s(e')$.

*Proof.* Let $i = \max(e \triangle e')$; $e <_{\text{colex}} e'$ means $i \in e'$, so $i \notin e$, and above
$i$ the sets agree. Let $m = p_i$.
(a) $m \in s(e')$ and $m \notin s(e)$: $m$ lies in $b_j$ only for $j = i$ (reduced), and
$i \in e'$, $i \notin e$.
(b) For a digit $x > m$: if $x \in b_j$ then $p_j \ge x > m = p_i$, so $j > i$; there $e$ and $e'$
agree, so $x$ lies in the same number of $b_j$ for $e$ and for $e'$: $x \in s(e) \iff x \in s(e')$.
So every digit in $s(e) \setminus s(e')$ is $< m$ (not $> m$ by (b), not $m$ by (a)), and
$m \in s(e') \setminus s(e)$: that is $s(e) <_{\text{colex}} s(e')$. $\square$

Since $k = \sum_{i \in e} 2^i$ orders the $e$ exactly as colex does (binary numbers; mathlib's
`geomSum_lt_geomSum_iff_toColex_lt_toColex`), $k \mapsto s(e)$ is strictly increasing in $k$
(Lean: `spanSum_lt_of_rank_lt`). It is onto $\mathbb{F}_{2^n}$ (a basis spans), so it is the
increasing enumeration: the element of rank $k$ is the $k$-th smallest. Linearity:
$s(e \triangle e') = s(e) \triangle s(e')$. $\square$

```lean
theorem spanSum_strictMono (b : ι → Finset α) (hne : ∀ i, (b i).Nonempty)
    (hp : StrictMono fun i => (b i).max' (hne i))
    (hred : ∀ i j, i ≠ j → (b i).max' (hne i) ∉ b j) :
    StrictMono fun e : Colex (Finset ι) => toColex (spanSum b (ofColex e))

theorem spanSum_lt_of_rank_lt (b : ℕ → Finset α) (hne : ∀ i, (b i).Nonempty)
    (hp : StrictMono fun i => (b i).max' (hne i))
    (hred : ∀ i j, i ≠ j → (b i).max' (hne i) ∉ b j) {e e' : Finset ℕ}
    (h : ∑ i ∈ e, 2 ^ i < ∑ i ∈ e', 2 ^ i) :
    toColex (spanSum b e) < toColex (spanSum b e')
```

## Not formalized

Step 1 (the digit structure of $\mathbf{On}_2$ below $\omega^{\omega^\omega}$) rests on DiMuro's
Lemma 1.1 and Theorem 3.1; connecting it to CombinatorialGames' Lean `Nimber` (which proves the
simplest extension theorems) is open. Comparison with $\mathbb{R}$, $\mathbb{C}$: there is no
analogue: an ordered field has no finite subfields, and here the order is not a field order
(it is compatible neither with multiplication nor with translation $x \mapsto x + c$); the lemma
says only that rank coordinates enumerate the subfield in increasing order.

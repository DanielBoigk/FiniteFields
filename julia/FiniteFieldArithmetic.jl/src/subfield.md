#implementation #characteristic-2 #theorem

> The smallest field of an element (backlog J9): `minimal_level(x)`, `minimal_degree(x)`
> $= [\mathbb{F}_2(x) : \mathbb{F}_2]$ and `reduce_field(x)`, `x` in that field. For nimbers all
> three come from the bit length of `x`, with no field arithmetic: `minimal_level` 0.1–0.6 ns,
> `reduce_field` 3.5–4.7 ns (type-unstable by design).

> Sources: R. Lidl, H. Niederreiter, *Finite Fields* (2nd ed., 1997), Theorem 2.6 (the subfields
> of $\mathbb{F}_{p^n}$ are the $\mathbb{F}_{p^d}$, $d \mid n$, one for each $d$) and Lemma 2.4 /
> Theorem 2.14 ($\mathbb{F}_{q^d}$ is the set of roots of $x^{q^d} - x$); J. H. Conway, *On Numbers
> and Games* (2nd ed., 2001), ch. 6 (the nimber fields, nested as integers).
> mathlib: [`FiniteField.isSplittingField_sub`](https://leanprover-community.github.io/mathlib4_docs/Mathlib/FieldTheory/Finite/GaloisField.html#FiniteField.isSplittingField_sub)
> (a field of $q$ elements is the splitting field of $X^q - X$); not formalized here.
> code: `src/subfield.jl`. tests: `test/subfield.jl`. benchmark: `bench/subfield.jl`.

## Theory

$\mathbb{F}_2(x)$ is the smallest subfield of $\mathbb{F}_{2^n}$ containing $x$. Its degree $d$ is
the degree of the minimal polynomial of $x$ and the smallest $d \ge 1$ with $x^{2^d} = x$, since
$\mathbb{F}_{2^d}$ is exactly the fixed set of the $d$-th power of the Frobenius. It divides $n$
(Lidl–Niederreiter, Thm. 2.6).

**Nimbers.** In $\mathbb{F}_{2^{2^K}}$ every $d \mid 2^K$ is a power of two, $d = 2^j$, and the
unique subfield of that degree is the nimber level $j$: the integers below $2^{2^j}$ (Conway,
ch. 6). So $\mathbb{F}_2(x)$ is the smallest level containing $x$ as an integer:

$$
\text{minimal\_level}(x) = \min\{\, j \ge 0 : x < 2^{2^j} \,\} = \lceil \log_2 \ell \rceil
\quad(\ell = \text{bit length of } x,\ \text{and } 0 \text{ for } \ell \le 1),
\qquad \text{minimal\_degree}(x) = 2^{\text{minimal\_level}(x)}.
$$

Compare $\mathbb{R} \subset \mathbb{C}$: there is no analogue of a "smallest field" read off a
representation; here the nesting of the levels as integers makes it a bit count.

## Implementation

- `minimal_level` computes $\ell$ with `leading_zeros` and $\lceil \log_2 \ell \rceil$ as
  `64 - leading_zeros(max(ℓ, 1) - 1)`, branch-free; `minimal_degree` shifts. Both are
  `@stable` and checked by JET.
- `reduce_field` returns the element in level `minimal_level(x)`, same bits. Its type depends on
  the value, so it is outside `@stable` (maintainer, 2026-11-20). A first version built the type
  `Nimber{j, nimber_storage(j)}` at runtime and cost 180–205 ns per element (dynamic dispatch);
  one branch per level, each with a concrete type, costs 3.5–4.7 ns. A caller still sees a Union
  of up to eight types, so inside a loop `minimal_level` (stable) is the better tool.
- Other representations (Conway, polynomial bases: J1a/J1b; $\overline{\mathbb{F}}_2$: J7) will
  need field arithmetic: the smallest $d \mid n$ with $x^{2^d} = x$, then the map down to
  $\mathbb{F}_{2^d}$ through the Conway compatibility. Not implemented.

## Tests

The oracle is the definition, independent of the bit-length shortcut: the smallest $d$ with
$x^{2^d} = x$, by squaring in the field. Checked on every element of GF(2) … GF(2^16) and on
random elements of GF(2^32) … GF(2^128), including random elements of each of their subfields and
the edges $2^{2^j} - 1$, $2^{2^j}$. `reduce_field` is checked to give the same element
(`==`, same bits), in the field of degree `minimal_degree(x)`, idempotently, with products and
quotients agreeing with those in the original field.

## Measurements

`bench/subfield.jl`, 1024 random elements, Julia 1.13.1, x86_64-linux (the development VM),
single run, ns per element:

| field | `minimal_level` | `reduce_field` |
|---|---|---|
| GF(256) | 0.30 | 4.7 |
| GF(2^16) | 0.28 | 3.8 |
| GF(2^32) | 0.18 | 3.7 |
| GF(2^64) | 0.14 | 4.3 |
| GF(2^128) | 0.58 | 3.5 |

Random elements almost always lie in the top level, so these mostly time one branch of
`reduce_field`.

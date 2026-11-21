#overview

> Handoff for the next session: where the work stands, what is next, and what is open. A session
> that starts reads it after `AGENTS.md`; a session that ends updates it. Replace stale text
> rather than appending; history belongs in git.

**Last updated:** 2026-11-21, on branch **`dev-danielb`** (the working branch; pushed).
`plan/agents-handoff` is merged into it and no longer used.

## Stage

Implementation of approved items. Approved and done so far: J0, the nimber fields GF(2) …
GF(2^128), J2a, J8 (maintainer, 2026-11-19: "please implement J8!"), J8c (2026-11-19: "I
kinda want all three"), J8d (2026-11-20: "go with your J8d proposal, fix the GF(2^16)
slowdown first"), J8e (2026-11-20: "Switch to the tables"), J2b, the default for `inv`
(2026-11-20: "Please implement the next item!", taken as approval of item 1 of the then "Next"
list), **J8b** (2026-11-20: "Please implement the next item!" again, taken the same way;
packing, cache blocking measured, bitsliced matrix product and vectors done, the rest open) and
**J8g** (2026-11-20: "yes. once you benchmarked it somewhat it should become the default by
matrix size!"), **FF-P2** (2026-11-20: "Go ahead with it!"), **J9** (2026-11-20: asked "What
is the next item? Please implement!", then chose J9 when asked which) and **J1a** (2026-11-21:
"Please implement the next item!", then chose J1a, Conway fields, when asked which). Ask before
starting a backlog section the maintainer has not approved.

## Done last

- **J1a, first part: polynomial-basis and Conway fields** (`src/polyfield.jl`, note
  `src/polyfield.md`, data `src/conway_table.jl`, tests `test/conway.jl`, benchmark
  `bench/polyfield.jl`). `PolyElement{N, T, F}` = $\mathbb{F}_2[x]/(x^N + F)$ for $N \le 128$;
  `poly_field(n, low)`; `conway_field(n)` for the 113 degrees $\le 128$ in Lübeck's table (not
  128); `GF256_AES`, `GF2_128_GHASH`; `embed` between Conway fields ($m \mid n$) and between
  Conway and nimber fields of degree $2^k \le 64$. Products: carry-less multiplication plus
  Barrett reduction (exactness proved in the note). Inverses: extended Euclid (Fermat at 128).
  A small refactor: `PolyBasis(F, α)` in `src/nimber_fast.jl` takes the root.
- **Measured** (`src/polyfield.md`): Conway products 2.5–3 ns up to degree 64, 9 ns above;
  *faster than the nimber defaults from GF(2^32) up* (2.5 vs 18.9 ns at 64; GHASH 8.4 vs 40.5 at
  128), slower below (2.9 vs 0.4 ns at GF(256)). Inverses slow (400 ns at 64, 3 µs GHASH).
  Nimber ↔ Conway 0.2–1.9 ns per element; Conway → Conway 335 ns – 1.5 µs.
- **Verified:** `Pkg.test()` passes (Julia 1.13.1, x86_64-linux), incl. Aqua and JET
  (`@test_opt` on products/inverses of six polynomial fields). Table: every entry irreducible
  and compatible with its divisors, primitive for $n \le 64$, the least such polynomial for
  $n \le 16$ (independent BigInt oracle). Products against the oracle on every Conway field,
  exhaustive up to $n = 10$; FIPS-197 AES vectors.
- **Not verified:** primitivity of the table entries above degree 64 (needs factorisations of
  $2^n - 1$); a GCM test vector for `GF2_128_GHASH` (bit order differs, documented); other Julia
  versions; CI; the jldoctest of `conway_field` (no Documenter run).

Before that: J9 (smallest field of a nimber), FF-P2 (`CLMulMatMul` for matrix × vector), J8b
(inlining `LazyMulBy`), J8g (`BitslicedMatMul` the default by size), the bitsliced matrix
product and vectors at GF(2^64)/GF(2^128).
Earlier: J8b cache blocking (measured, not adopted), J8b packing by bit-matrix transposes, J2b,
J8e and J8d, J8c, J8, J2a, nimber GF(2) … GF(2^128), ADR-0004, `AGENTS.md`, CI, the move to
MathStruct. Details in git and the Backlog's *Done*.

## Next

Within the approved J1a (Backlog, "J1a — Representations, remaining"): byte-table Conway →
Conway embeddings, table strategies for small Conway fields, bulk operations for `PolyElement`
can go ahead. Two parts need a decision first: **pseudo-Conway polynomials** for $n = 128$ and the
other missing degrees (not canonical: which construction?), and **isomorphisms between AES/GHASH
and Conway fields** (one of $n$, or all). Within J8b only ARM `PMULL` and GPU kernels remain
(not doable in this container).

Not approved (need the maintainer's go): **J9b** (smallest field for Conway fields; now possible,
J1a exists); **J3** (inversion strategies; the polynomial-basis inverse needs it most: 400 ns at
degree 64); **J8f**; a cache-blocked transpose in `CLMulMatMul`; faster small bitsliced products;
**FF-P1**; **J1b** (arbitrary $n$, above 128); linear algebra (V5). Seen earlier, not measured:
in the scalar GF(2^16) `LazyMulBy` loop the size of the global `NIM256_MUL` table is reloaded
every element.

## Decisions in force (maintainer, 2026-11-19)

- Purpose: theory, cryptanalysis, and learners built on finite fields; no cryptographic
  implementations, constant time not a goal.
- Characteristic 2 first; byte-sized fields first ($\mathbb{F}_2$, $\mathbb{F}_{2^8}$ … $\mathbb{F}_{2^{128}}$),
  then arbitrary $n$ (J1b) and $\overline{\mathbb{F}}_2$ (J7).
- Representations: [[ADR-0004 Nimbers as the default representation of GF(2^n)]].
- Strategies: tables and computed variants side by side, the default chosen by benchmark.
- Printing: hexadecimal always; `show` adds the field when it does not fill its storage.
- Layout `julia/FiniteFieldArithmetic.jl/`; `lake-manifest.json` and `flake.lock` committed; CI
  path-filtered plus the `full-ci` branch; mathlib upstreaming case by case; milestones deferred.

## Environment notes

- Cloud sessions: use the personal environment **FiniteFields** (allows `release.lean-lang.org`,
  `releases.lean-lang.org`, `cache.mathlib.org`, `pkg.julialang.org`). Start sessions with
  `MathStruct/FiniteFields` as the source; this session still pushes through GitHub's redirect from
  the old `DanielBoigk/FiniteFields` URL.
- No system Julia in the container: `nix develop .#julia -c julia ...` fetches Julia 1.13.1.
  The bench environment uses `[sources]` (Julia ≥ 1.11).
- On 2026-11-20 the session's network refused `pkg.julialang.org` (registry update failed); tests
  and benchmarks still ran from the packages already in `~/.julia`. A new dependency would need
  the **FiniteFields** environment.

## Open questions for the maintainer

0. **J1a decisions taken provisionally (2026-11-21), please confirm or change:**
   (a) the nimber ↔ Conway isomorphism sends $\alpha_{2^k}$ to the *smaller integer* of the two
   compatible roots ($r_3 =$ `0x61` in GF(256); table in `src/polyfield.md`); (b) no implicit
   promotion or `==` between fields: maps only through `embed`; (c) names `PolyElement`,
   `poly_field`, `conway_field`, `embed`, `GF256_AES`, `GF2_128_GHASH`; (d) the Conway data is
   Lübeck's, taken from the GPL-3 `sagemath/conway-polynomials` mirror (provenance and sha256 in
   `src/conway_table.jl`); the polynomials are mathematical facts, but say if you want the
   licensing checked or a different source; (e) Conway/polynomial-basis products are 4–7× faster
   than the nimber defaults from GF(2^32) up (measurement in `src/polyfield.md`): does that
   change anything in ADR-0004 for bulk work?
1. **Learners on finite fields:** when you are ready, write the ideas down (a research-question
   note in `vault/Research Questions/`), so the arithmetic API can be checked against them.
2. **GitHub Pages:** check that the repository's Pages setting serves the `gh-pages` branch.
3. **`main`:** merged on 2026-11-19 at the maintainer's request (fast-forward to `40cd1c2`, J2a).
   That push should run `julia.yml`, `lean.yml` and `pages.yml`; their results could not be
   read from the session (no GitHub API access), so check the Actions tab and the site.
4. **Milestones** (deferred).

## Known issues

- `lake build` warns that `lean/FiniteFields.lean` should start with `module` (mathlib uses the
  module system). Adopt it with the Lean layout (L0).
- Wikilinks use two styles (`[[Decisions/…]]` and `[[Project/Decisions/…]]`); move to by-name links
  with V0.

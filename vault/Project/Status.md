#overview

> Handoff for the next session: where the work stands, what is next, and what is open. A session
> that starts reads it after `AGENTS.md`; a session that ends updates it. Replace stale text
> rather than appending; history belongs in git.

**Last updated:** 2026-11-20, on branch **`dev-danielb`** (the working branch; pushed).
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
matrix size!"), **FF-P2** (2026-11-20: "Go ahead with it!") and **J9** (2026-11-20: asked "What
is the next item? Please implement!", then chose J9 when asked which). Ask before starting a backlog
section the maintainer has not approved.

## Done last

- **J9, the smallest field of an element (nimbers)**: `minimal_level`, `minimal_degree`,
  `reduce_field` (exported) in `src/subfield.jl`, note `src/subfield.md`, benchmark
  `bench/subfield.jl`. All from the bit length, since every subfield of
  $\mathbb{F}_{2^{2^K}}$ is a nimber level. `reduce_field` is type-unstable (allowed by the
  maintainer); one branch per level made it 3.5–4.7 ns per element instead of 180–205 ns when
  it built `Nimber{j, S}` at runtime. The other representations are J9b (after J1a/J1b/J7).
- **Verified:** `Pkg.test()` passes (947,134 checks, Julia 1.13.1, x86_64-linux), including the
  Frobenius oracle on every element of GF(2) … GF(2^16), JET `@test_opt` on `minimal_level` and
  `minimal_degree`; `bench/subfield.jl` (single run).
- **Not verified:** other Julia versions; CI; the jldoctest in `minimal_degree`'s docstring
  (no Documenter run; checked by hand: 0x1b has 5 bits, level 3, degree 8).

Before that: FF-P2, `CLMulMatMul` for matrix × vector (a column-wise path chosen by the size of
$A$, `src/bulk.md`, "Matrix × vector"); J8b, why inlining both branches of `LazyMulBy`'s
`copyto!` slowed the loop; J8g (`BitslicedMatMul` the default by size); the bitsliced matrix
product and vectors at GF(2^64)/GF(2^128).
Earlier: J8b cache blocking (measured, not adopted), J8b packing by bit-matrix transposes, J2b,
J8e and J8d, J8c, J8, J2a, nimber GF(2) … GF(2^128), ADR-0004, `AGENTS.md`, CI, the move to
MathStruct. Details in git and the Backlog's *Done*.

## Next

Within the approved J8b, only ARM `PMULL` and GPU kernels remain: they need an ARM machine and a
GPU and are not doable in this container. Everything else needs the maintainer's go.

Not approved (need the maintainer's go): **J8f** (fewer carry-less products in `CLMulMatMul`);
a **cache-blocked transpose** in `CLMulMatMul`'s transposed path (might make the column path
unnecessary and help small $q \ge 3$; an expectation, `src/bulk.md`); faster small bitsliced
products (preallocated temporaries, a `BitslicedMatrix` type keeping the planes; would lower the
thresholds); **FF-P1** (`CLMul` at GF(2^64) element by element); **J1a/J1b** (Conway
representation, arbitrary $n$); **J3** (Itoh–Tsujii, extended Euclid as `inv` strategies);
linear algebra on the matrix products (V5); **J9b** (smallest field for Conway and polynomial
representations, after J1a/J1b/J7). Seen earlier, not measured: in the scalar GF(2^16)
`LazyMulBy` loop the size of the global `NIM256_MUL` table is reloaded every element.

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

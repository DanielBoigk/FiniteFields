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
list), and **J8b** (2026-11-20: "Please implement the next item!" again, taken the same way).
Within J8b: packing done, cache blocking measured and not adopted, and now (the fourth "next
item" of 2026-11-20: "Implement the next item!") the bitsliced matrix product and bitsliced
vectors at every level. J8b's remaining sub-items are listed under "Next".
Ask before starting a backlog section the maintainer has not approved.

## Done last

- **J8b, bitsliced matrix product — `BitslicedMatMul`** (`src/bitsliced_matmul.jl`, note
  `src/bitsliced_matmul.md`): a matrix over GF(2^n) as n bit matrices; the nimber recursion on
  whole bit matrices ($3^K$ GF(2) products, correct because every step is bilinear or linear),
  each GF(2) product by the Four Russians (eight byte tables per 64 inner indices, 256-bit
  accumulators). `bench/matmul_sizes.jl compare`, conversions included, from 1024 × 1024:
  90–170× faster than the default at GF(2) … GF(16), 48× at GF(256), 7–13× at GF(2^16),
  2.8–3.6× at GF(2^32), 2.1–2.2× at GF(2^64), 2.9× at GF(2^128); slower on small matrices
  (crossovers 64 … 512 by field). **Not made a default**: proposed as **J8g** (choose by size).
- **J8b, `BitslicedVector` at GF(2^64) and GF(2^128)** (`src/bulk.md`): the top one or two
  recursion levels through a scratch buffer, the level-5 circuit unrolled: 2.6 and 7.4 ns per
  element against 18.4 and 41.8 for the scalar product (scratch run); GF(2^128) packs its
  64-bit halves separately (6.0 ns per element).
- **Verified:** `Pkg.test()` passes, 733,010 checks, Julia 1.13.1, x86_64-linux, Aqua and JET
  included (new: every matrix-product test with `BitslicedMatMul`, block-boundary shapes,
  linearity, a 300 × 1100 × 70 product at every level, bitsliced vectors at every level with an
  allocation bound). `bench/matmul_sizes.jl compare` once (the table in the note).
- **Not verified:** CI (not pushed to `full-ci`); Julia 1.10; ARM; `bench/benchmarks.jl` with the
  new rows was not run as a whole; the GF(2^64)/GF(2^128) vector numbers come from a scratch
  script, not from the committed benchmark.

Earlier: J8b cache blocking (measured, not adopted), J8b packing by bit-matrix transposes, J2b,
J8e and J8d, J8c, J8, J2a, nimber GF(2) … GF(2^128), ADR-0004, `AGENTS.md`, CI, the move to
MathStruct. Details in git and the Backlog's *Done*.

## Next

Within the approved J8b (remaining sub-items):

1. Why inlining both branches of `LazyMulBy`'s `copyto!` slowed the loop (J8d, not inspected).
2. ARM `PMULL` and GPU kernels: need an ARM machine and a GPU; not doable in this container.

Not approved (need the maintainer's go): **J8g** (`BitslicedMatMul` as the default from a
measured size per field; the measurements suggest it is the largest available speed-up for
matrix products); **J8f** (fewer carry-less products in `CLMulMatMul`, for the small sizes where
it stays the default); **FF-P1** (`Problems.md`, `CLMul` at GF(2^64) element by element);
**J1a/J1b** (Conway representation, arbitrary $n$); **J3** (Itoh–Tsujii, extended Euclid as
`inv` strategies); linear algebra on the matrix products (V5).

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

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
matrix size!"). Ask before starting a backlog section the maintainer has not approved.

## Done last

- **J8g, `BitslicedMatMul` as the default by size.** `default_matmul_strategy(F)` is
  `BySizeMatMul(TableMatMul() or CLMulMatMul(), BitslicedMatMul(), bitsliced_threshold(F))`:
  bitsliced when every dimension of the product reaches the threshold, 48 (GF(2), GF(4)), 64
  (GF(16)), 96 (GF(256)), 192 (GF(2^16)), 256 (GF(2^32)), 512 (GF(2^64)), 384 (GF(2^128)).
  Measured with the new `bench/matmul_threshold.jl` (square products on a fine grid, then
  rectangular shapes: with the smallest dimension at the threshold the bitsliced product won in
  every position, 1.4–20×; with a dimension of 1 or 2 it loses and is not chosen). Details in
  `src/bitsliced_matmul.md`, "Default by size".
- **Found on the way:** FF-P2 in `Problems.md`, `CLMulMatMul` is slow for matrix × vector (9–70
  ns per multiply-add at GF(2^32) … GF(2^128) for 1024 × 1024 times a vector). Pre-existing, not
  fixed.
- Before that, in the same session: `BitslicedMatMul` itself and `BitslicedVector` at GF(2^64) and
  GF(2^128) (J8b; `src/bitsliced_matmul.md`, `src/bulk.md`).
- **Verified:** `Pkg.test()` passes (733,299 checks), Julia 1.13.1, x86_64-linux, Aqua and
  JET included; new tests for the choice by size, `*` and `mul!` on both sides of every
  threshold. `bench/matmul_threshold.jl` once.
- **Not verified:** CI (not pushed to `full-ci`); a second threshold run (the VM varies up to
  about 40%, so each threshold is uncertain by one grid step; at GF(2^32) the margin is 1.05×);
  other CPUs (thresholds measured on one Intel cloud VM); Julia 1.10; ARM.

Earlier: J8b cache blocking (measured, not adopted), J8b packing by bit-matrix transposes, J2b,
J8e and J8d, J8c, J8, J2a, nimber GF(2) … GF(2^128), ADR-0004, `AGENTS.md`, CI, the move to
MathStruct. Details in git and the Backlog's *Done*.

## Next

Within the approved J8b (remaining sub-items):

1. Why inlining both branches of `LazyMulBy`'s `copyto!` slowed the loop (J8d, not inspected).
2. ARM `PMULL` and GPU kernels: need an ARM machine and a GPU; not doable in this container.

Not approved (need the maintainer's go): **FF-P2** (a matrix × vector path for `CLMulMatMul`);
**J8f** (fewer carry-less products in `CLMulMatMul`, which is now the default only below the
bitsliced thresholds); faster small bitsliced products (preallocated temporaries, a
`BitslicedMatrix` type keeping the planes; would lower the thresholds); **FF-P1** (`CLMul` at
GF(2^64) element by element); **J1a/J1b** (Conway representation, arbitrary $n$); **J3**
(Itoh–Tsujii, extended Euclid as `inv` strategies); linear algebra on the matrix products (V5).

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

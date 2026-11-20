#overview

> Handoff for the next session: where the work stands, what is next, and what is open. A session
> that starts reads it after `AGENTS.md`; a session that ends updates it. Replace stale text
> rather than appending; history belongs in git.

**Last updated:** 2026-11-20, on branch **`dev-danielb`** (the working branch; pushed).
`plan/agents-handoff` is merged into it and no longer used.

## Stage

Implementation of approved items. Approved and done so far: J0, the nimber fields GF(2) …
GF(2^128), J2a, J8 (maintainer, 2026-11-19: "please implement J8!") and J8c (2026-11-19: "I
kinda want all three"). Ask before starting a
backlog section the maintainer has not approved.

## Done last

- **J8c — strategies for `c .* a`** (`src/bulk.jl`, note `src/bulk.md`, "Constant times array"),
  tests and benchmark first (`test/constmul.jl`, `bench/constmul.jl`):
  - `c .* a`, `a .* c` and `c * a` go through `default_constmul_strategy(F)`: `ScalarConstMul()`
    (C, Base's behaviour, **still the default everywhere**), `EagerMulBy()` (A), `LazyMulBy(t)`
    (B: a lazy broadcast of a callable that builds `MulBy` tables only for ≥ `t` elements).
    Switch per field by redefining `default_constmul_strategy`.
  - Measured break-even of B (now `mulby_threshold`): 768 / 192 / 256 / 512 elements at
    GF(2^16) / GF(2^32) / GF(2^64) / GF(2^128); higher than J8's estimate of 40–240. On long arrays
    B is 2.4–9 times faster than C from GF(2^16) up, alone and fused.
  - At GF(256) B gains about 2× in `c .* a` and `out .= c .* a`, nothing in `out .= c .* a .+ b`;
    there A (eager) wins the fused case (0.25–0.30 vs 0.48 ns), not profiled.
  - Weakness found: B below its threshold is 19–27% slower than C at GF(2^16) in fused
    expressions (per-element branch, presumably); elsewhere equal.
  - `ConstMul` uses its own 256-byte table at GF(256) too: 0.20 vs 0.36 ns for indexing the
    shared product table.
- **Verified:** `Pkg.test()` passes, 733,795 checks, Julia 1.13.1, x86_64-linux, Aqua and JET
  included. Benchmarks: one full run on the noisy cloud VM (outliers visible in the raw output;
  thresholds read off the first length from which B won every expression).
- **Not verified:** Julia 1.10; ARM; CI (not pushed to `full-ci`); GF(2) … GF(16) thresholds
  (taken from GF(256)); the cause of the GF(256) and GF(2^16) effects (not profiled).

Earlier: J8 (bulk operations: `rand!`, `MulBy`, matrix products, `BitslicedVector`), J2a, nimber
GF(2) … GF(2^128), ADR-0004, `AGENTS.md`, CI, the move to MathStruct. Details in git.

## Next (proposals, need the maintainer's go)

1. **Defaults for `inv`:** `LogTable16` inverts about 2.3 times faster at GF(2^16) and 1.9 times at
   GF(2^32), but its products have about twice the latency of `MulTable`. A separate default
   for `inv` (`default_inv_strategy`) would take both. Decide whether that API is wanted.
2. **J8d — default for `c .* a`: decide.** Proposal (`src/bulk.md`): `LazyMulBy(F)` from GF(2^16)
   up, `ScalarConstMul` at GF(256) and below, drop `EagerMulBy`; fix B's 20% below-threshold
   cost at GF(2^16) first. Nothing changes until the maintainer says so.
3. **J8b** (Backlog): faster bitsliced packing, cache blocking for large matrices, bitsliced
   matrix products, GPU kernels.
4. **FF-P1** (`Problems.md`): `CLMul` at GF(2^64) element by element; investigate on a quiet
   machine. (Matrix products at GF(2^64) no longer depend on it.)
5. **ARM `PMULL`** for `CLMul` and `CLMulMatMul` (needs an ARM machine to test).
6. **J1a/J1b** — Conway representation and arbitrary $n$ (e.g. GF(2^7)); J3 (inversion), and
   linear algebra on top of the matrix products (rank, solving; V5).

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

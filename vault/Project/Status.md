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
slowdown first") and J8e (2026-11-20: "Switch to the tables"). Ask before starting a
backlog section the maintainer has not approved.

## Done last

- **J8e — `c .* a` uses `MulBy` tables at GF(256) and below too:** `LazyMulBy(F)` is now the
  default for every field; threshold 96 at GF(2), GF(4), GF(16), GF(256), measured with
  `bench/constmul.jl small` (new option). Long arrays: 1.6–2.3× faster in `out .= c .* a`,
  1.3–1.4× in `out .= c .* a .+ b`. `ScalarConstMul` stays available but is no default.
- **J8d — `c .* a` uses `MulBy` tables by default from GF(2^16)** (`src/bulk.jl`, note
  `src/bulk.md`, "Constant times array"), tests first (`test/constmul.jl`):
  - `default_constmul_strategy`: `LazyMulBy(F)` from GF(2^16) up (thresholds 768 / 192 / 256 /
    512 at GF(2^16) … GF(2^128), unchanged after three runs); GF(256) and below followed in J8e.
    `EagerMulBy` and the J8c callable `ConstMul` are removed.
  - The GF(2^16) slowdown is fixed: `LazyMulBy` is now a broadcast style (`ConstMulStyle`) whose
    marker spreads to the fused expression; `copyto!`/`copy` decide once by the result's length
    and rewrite the tree, each branch behind a `@noinline` barrier. Below the threshold it equals
    Base's loop from 1024 elements (1.50 vs 1.50 ns), a few ns per call more below. Two simpler
    designs (per-element flag, `Union` return) measured 1.80 and 1.74 ns; recorded in the note.
  - From GF(2^16) up `c .* a` is now 2.4–9 times faster on long arrays, fused included.
  - The maintainer's question about an allocation-free eager `EagerMul!`: measured
    (`out .= m.(a); out .+= b` with a prepared `MulBy`), pays only at GF(256) in cache, loses
    beyond the caches from GF(2^32) up; no new API, the answer is in the note.
- **Verified:** `Pkg.test()` passes, 732,509 checks, Julia 1.13.1, x86_64-linux, Aqua and JET
  included. Benchmarks: three full runs of `bench/constmul.jl` and one with `small`, on the noisy
  VM.
- **Not verified:** Julia 1.10; ARM; CI (not pushed to `full-ci`); the GF(2) … GF(16) thresholds
  rest on one run; why inlining both branches slowed the loop (not inspected); interaction of `ConstMulStyle` with
  other packages' broadcast styles (GPU arrays, StaticArrays): it wraps their style and hands the
  rewritten tree back to it, tested only with Base arrays.

Earlier: J8c (strategies and measurements for `c .* a`), J8 (bulk operations), J2a, nimber
GF(2) … GF(2^128), ADR-0004, `AGENTS.md`, CI, the move to MathStruct. Details in git.

## Next (proposals, need the maintainer's go)

1. **Defaults for `inv`:** `LogTable16` inverts about 2.3 times faster at GF(2^16) and 1.9 times at
   GF(2^32), but its products have about twice the latency of `MulTable`. A separate default
   for `inv` (`default_inv_strategy`) would take both. Decide whether that API is wanted.
2. **J8b** (Backlog): faster bitsliced packing, cache blocking for large matrices, bitsliced
   matrix products, GPU kernels.
3. **FF-P1** (`Problems.md`): `CLMul` at GF(2^64) element by element; investigate on a quiet
   machine. (Matrix products at GF(2^64) no longer depend on it.)
4. **ARM `PMULL`** for `CLMul` and `CLMulMatMul` (needs an ARM machine to test).
5. **J1a/J1b** — Conway representation and arbitrary $n$ (e.g. GF(2^7)); J3 (inversion), and
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

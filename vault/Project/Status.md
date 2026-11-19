#overview

> Handoff for the next session: where the work stands, what is next, and what is open. A session
> that starts reads it after `AGENTS.md`; a session that ends updates it. Replace stale text
> rather than appending; history belongs in git.

**Last updated:** 2026-11-19, on branch **`dev-danielb`** (the working branch; pushed).
`plan/agents-handoff` is merged into it and no longer used.

## Stage

Implementation of approved items. Approved and done so far: J0, the nimber fields GF(2) …
GF(2^128), and J2a (maintainer, 2026-11-19: "Implement J2a"). Ask before starting a backlog
section the maintainer has not approved.

## Done last

- **J2a — faster large nimber levels** (`src/nimber_fast.jl`, note `src/nimber_fast.md`), tests and
  benchmarks first. Two new strategies:
  - `LogTable16`: the recursion stops at GF(2^16), with log/exp tables (384 KiB).
  - `CLMul`: a polynomial basis GF(2)[x]/(f) for sparse irreducible f (GHASH polynomial at 128
    bits), reached through a change of basis (byte tables). The root of f in the nimber field is
    found at precompile time (Cantor–Zassenhaus). The product uses x86 `PCLMULQDQ` through an LLVM
    intrinsic and is used at GF(2^64) and GF(2^128). Without the instruction (ARM, very old x86) it
    computes as `MulTable`.
- **Default changed:** `default_strategy(GF2_128) = CLMul()`. In throughput it measured 46–53 ns
  per product against 77–80 ns for `MulTable` (latency 51–55 against 72–74), in every run. The
  other levels keep `MulTable`.
- **Verified:** `Pkg.test()` passes, 728,712 checks, Julia 1.13.1, x86_64-linux with `PCLMULQDQ`.
  The tests also exercise the path without the instruction by switching the flag off. Aqua and
  JET pass.
- **Not verified:** Julia 1.10 (`min` in CI); ARM; an x86 CPU without `PCLMULQDQ`; CI as a whole
  (not pushed to `full-ci` from this session).
- **Benchmark noise:** this VM varies by up to about 40% between runs; `src/nimber_fast.md` gives
  ranges over several runs.

Earlier: nimber GF(2) … GF(2^128) (`Nimber{0..7}`), ADR-0004, `AGENTS.md`, CI, the move to
MathStruct. Details in git.

## Next (proposals, need the maintainer's go)

1. **Defaults for `inv`:** `LogTable16` inverts about 2.3 times faster at GF(2^16) and 1.9 times at
   GF(2^32), but its products have about twice the latency of `MulTable`. A separate default
   for `inv` (`default_inv_strategy`) would take both. Decide whether that API is wanted.
2. **FF-P1** (`Problems.md`): `CLMul` at GF(2^64) is about 1.8 times faster than `MulTable` on bare
   bit patterns but not consistently so as `GF2_64` elements. Investigate on a quiet machine, then
   reconsider the GF(2^64) default.
3. **J8 — bulk multiplication.** With data kept in the polynomial basis, a GF(2^128) product
   costs 4.2 ns instead of about 46. Converting once per array is the obvious bulk path, and it
   is also the route to SIMD/GPU.
4. **ARM `PMULL`** for `CLMul` (needs an ARM machine to test).
5. **J1a/J1b** — Conway representation and arbitrary $n$ (e.g. GF(2^7)).

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

## Open questions for the maintainer

1. **Learners on finite fields:** when you are ready, write the ideas down (a research-question
   note in `vault/Research Questions/`), so the arithmetic API can be checked against them.
2. **GitHub Pages:** check that the repository's Pages setting serves the `gh-pages` branch.
3. **`main`:** `dev-danielb` holds everything since 2026-11-19; merge it into `main` when you
   want the site deployed (`pages.yml` deploys from `main`).
4. **Milestones** (deferred).

## Known issues

- `lake build` warns that `lean/FiniteFields.lean` should start with `module` (mathlib uses the
  module system). Adopt it with the Lean layout (L0).
- Wikilinks use two styles (`[[Decisions/…]]` and `[[Project/Decisions/…]]`); move to by-name links
  with V0.

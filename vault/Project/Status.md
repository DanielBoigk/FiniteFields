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
slowdown first"), J8e (2026-11-20: "Switch to the tables") and J2b, the default for `inv`
(2026-11-20: "Please implement the next item!", taken as approval of item 1 of the previous
"Next" list). Ask before starting a backlog section the maintainer has not approved.

## Done last

- **J2b — a separate default for `inv`:** `default_inv_strategy(F)` (`src/nimber.jl`,
  `src/nimber_fast.jl`) is behind `inv` and so `/`; it equals `default_strategy(F)` except
  `LogTable16()` at GF(2^16) and GF(2^32). `*` is unchanged. Tests first (`test/inv_default.jl`,
  JET on `inv` and `/`); the suite gained `inv latency` chains and `default (inv)` entries.
  Measurements and reasoning in `src/nimber_fast.md`, "Default for inv":
  - GF(2^16): `inv` 2.7× the throughput (0.85 vs 2.3–2.5 ns) and ~30% less latency.
  - GF(2^32): **a trade-off for the maintainer to confirm.** 1.4–1.9× the throughput (8–11 vs
    15 ns) but 20–50% more latency in a chain of dependent inversions (39–51 vs 31–33 ns).
    Throughput was chosen (the bulk case); say if you prefer `MulTable` there.
  - GF(2^64), GF(2^128): unchanged, `LogTable16` is not faster there.
- **Verified:** `Pkg.test()` passes, 732,584 checks, Julia 1.13.1, x86_64-linux, Aqua and JET
  included. Benchmarks: the `inv` groups, three runs in one process, on the noisy VM.
- **Not verified:** Julia 1.10; ARM; CI (not pushed to `full-ci`); the full benchmark suite was
  not rerun (only the `inv` groups); why `LogTable16`'s inversion has more latency at GF(2^32)
  (not profiled).

Earlier: J8e and J8d (`MulBy` tables behind `c .* a` at every field), J8c, J8 (bulk
operations), J2a, nimber GF(2) … GF(2^128), ADR-0004, `AGENTS.md`, CI, the move to MathStruct.
Details in git and the Backlog's *Done*.

## Next (proposals, need the maintainer's go)

1. **J8b** (Backlog): faster bitsliced packing, cache blocking for large matrices, bitsliced
   matrix products, GPU kernels.
2. **FF-P1** (`Problems.md`): `CLMul` at GF(2^64) element by element; investigate on a quiet
   machine. (Matrix products at GF(2^64) no longer depend on it.)
3. **ARM `PMULL`** for `CLMul` and `CLMulMatMul` (needs an ARM machine to test).
4. **J1a/J1b** — Conway representation and arbitrary $n$ (e.g. GF(2^7)); J3 (inversion:
   Itoh–Tsujii and extended Euclid as further `inv` strategies, now that `inv` has its own
   default), and linear algebra on top of the matrix products (rank, solving; V5).

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

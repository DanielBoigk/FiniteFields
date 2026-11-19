#overview

> Handoff for the next session: where the work stands, what is next, and what is open. A session
> that starts reads it after `AGENTS.md`; a session that ends updates it. Replace stale text
> rather than appending; history belongs in git.

**Last updated:** 2026-11-19, on branch `plan/agents-handoff` (from `main` at `a95d81f`; pushed;
not merged).

## Stage

**Implementation has started** (maintainer's go, 2026-11-19: "start with J0 and nimber GF(2^8)").
The planning-phase restriction in `AGENTS.md` is lifted for the backlog items the maintainer has
approved; ask before starting a new section of the backlog.

## Done last

- **J0 — `julia/FiniteFieldArithmetic.jl` skeleton**, test first (the tests were committed failing,
  then the implementation): `test/` (Aqua.jl, JET.jl `test_package` and `@test_opt`), DispatchDoctor
  `@stable` on all arithmetic, `bench/` (BenchmarkTools `SUITE`, `run.jl`, `profile_nimber.jl`),
  implementation notes `src/binaryfield.md`, `src/nimber.md`. CI: `julia.yml` now tests Julia `1`
  and `min` (1.10, the `[compat]` floor).
- **Nimber fields GF(2) ⊂ GF(4) ⊂ GF(16) ⊂ GF(256)** (`Nimber{K, UInt8}`, aliases `GF2` … `GF256`),
  per [[ADR-0004 Nimbers as the default representation of GF(2^n)]]. Strategies `Recursive`,
  `MulTable`, `LogTable`; default `MulTable` (0.44–0.49 ns per product in bulk, 2.1 ns latency).
  Hex printing; the REPL names the field when it does not fill its byte.
- **Verified:** `Pkg.test()` passes, 265,388 checks, on Julia 1.13.1 (x86_64-linux): the mex
  definition below 16, the Fermat-power rules, exhaustive field axioms on GF(256), exhaustive
  agreement of the strategies, Frobenius and the cyclic group, subfield nesting, printing, Aqua,
  JET. Benchmarks and profile recorded in `src/nimber.md`.
- **Not verified:** Julia 1.10 (`min`, only in CI, which has not run); macOS/ARM; GPU use.
- A bug found while checking a note: `GF256(7) == 7` threw; equality with other numbers is now
  `false`, with a regression test.

Earlier the same day: `AGENTS.md`, the flake fix (lazy Julia/elan), `flake.lock`,
`lean/lake-manifest.json` (`lake build` passes), path-filtered CI with a `full-ci` branch, the
move to `MathStruct/FiniteFields` (`baseUrl` `mathstruct.org/FiniteFields`). Details in git.

## Next

1. **J1 — larger nimber levels**: GF(2^16) (`UInt16`) … GF(2^128) (`UInt128`), recursing down to
   the GF(256) tables (a full table stops being possible at 16 bits). Tests: the same exhaustive
   checks are impossible above 16 bits, so use the subfield nesting (products of GF(256) values
   agree), random field-axiom checks, and the Fermat-power rules at every level.
2. **J8 — bulk multiplication**: `Recursive` is scalar (179 instructions, no 8-bit SIMD shifts on
   x86); bitslicing is the candidate for fast table-free and GPU arithmetic.
3. Then J1a/J1b (Conway representation, arbitrary $n$).

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
3. **Merging:** `plan/agents-handoff` holds everything since 2026-11-19; open a pull request to
   `main`?
4. **Milestones** (deferred).

## Known issues

- `lake build` warns that `lean/FiniteFields.lean` should start with `module` (mathlib uses the
  module system). Adopt it with the Lean layout (L0).
- Wikilinks use two styles (`[[Decisions/…]]` and `[[Project/Decisions/…]]`); move to by-name links
  with V0.

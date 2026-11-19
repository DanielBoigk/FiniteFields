#overview

> Handoff for the next session: where the work stands, what is next, and what is open. A session
> that starts reads it after `AGENTS.md`; a session that ends updates it. Replace stale text
> rather than appending; history belongs in git.

**Last updated:** 2026-11-19, on branch **`dev-danielb`** (the working branch; pushed).
`plan/agents-handoff` is merged into it and no longer used.

## Stage

Implementation of approved items. Approved and done so far: J0 and the nimber fields GF(2) …
GF(2^128). Ask before starting a backlog section the maintainer has not approved.

## Done last

- **Nimber fields GF(2^16), GF(2^32), GF(2^64), GF(2^128)** (`Nimber{4..7}` in `UInt16` … `UInt128`,
  aliases `GF2_16` … `GF2_128`), tests first. All strategies now share one recursion and differ in
  their base case (`Recursive`: GF(2); `MulTable`, `LogTable`: the GF(256) tables); halves are
  narrowed to the type of the level below. Powers reduce the exponent in the storage type (also
  `UInt128`, `BigInt` exponents).
- **Measured** (`src/nimber.md`): `MulTable` is fastest at every level and the default; ns per
  product in bulk: GF(256) 0.49, GF(2^16) 1.7, GF(2^32) 7.0, GF(2^64) 24.9, GF(2^128) 73. Profile:
  no allocations, no hotspot; half the GF(2^128) time is the 81+ base lookups.
- **Verified:** `Pkg.test()` passes, 463,506 checks, Julia 1.13.1, x86_64-linux (exhaustive on
  GF(256) and on GF(2^16) inverses; structural oracles above).
- **Not verified:** Julia 1.10 and CI as a whole (the `full-ci` run of the previous commit could
  not be read from this session; check the Actions tab); macOS/ARM; GPU.

Earlier: J0 and nimber GF(256) (ADR-0004); `AGENTS.md`; flake fix; lock files; path-filtered CI
with a `full-ci` branch; move to `MathStruct/FiniteFields`. Details in git.

## Next (proposals, need the maintainer's go)

1. **J2 for nimbers, speed of the large levels:** GF(2^128) costs 73 ns per product. Options:
   (a) GF(2^16) log/exp tables as the recursion base (one level fewer, ~256 KiB of tables);
   (b) carry-less multiplication in a polynomial basis behind a change of basis (a bit matrix),
   attractive in bulk. Measure both against the current default.
2. **J8 — bulk multiplication** (bitslicing; also the route to GPU).
3. **J1a/J1b** — Conway representation and arbitrary $n$ (e.g. GF(2^7)).

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

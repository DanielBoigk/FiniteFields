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

- **J8b, why inlining both branches of `LazyMulBy`'s `copyto!` slowed the loop** (taken as
  approved: item 1 of "Next" within the approved J8b, maintainer 2026-11-20 "Please implement the
  next item!"). Cause: LLVM's loop unswitching. Base's broadcast checks per element and argument
  whether the argument is extruded (a `select`); LLVM normally hoists both checks out of the loop,
  but its unswitching cost multiplier grows with the number of loops in a method, so with the
  table loop in the same method the scalar loop kept one check per element (GF(2^16), 4096
  elements fused: 1.384 ns against 1.254 for Base). With
  `JULIA_LLVM_ARGS=--enable-unswitch-cost-multiplier=false` the inlined variant is unswitched
  fully and runs at 1.255, which confirms it. GF(256) is not affected. The function barrier stays;
  no library code changed except a comment. Script `bench/constmul_inline.jl`, details in
  `src/bulk.md`, "Why inlining both branches slowed the loop".
- **Seen on the way, not measured:** in the scalar GF(2^16) loop the size of the global
  `NIM256_MUL` table is reloaded every element (possible aliasing with the store into `out`).
  Could be worth a look; not recorded as a problem because its cost is unknown.
- **Verified:** `bench/constmul_inline.jl` with and without the LLVM flag (GF(2^16)), once at
  GF(256); `Pkg.test()` passes (733,299 checks, Julia 1.13.1, x86_64-linux).
- **Not verified:** other fields than GF(2^16) and GF(256); other Julia/LLVM versions; CI.

Before that: J8g (`BitslicedMatMul` the default by size, thresholds 48 … 512 per field,
`bench/matmul_threshold.jl`), FF-P2 found (`CLMulMatMul` slow for matrix × vector, not fixed),
the bitsliced matrix product and vectors at GF(2^64)/GF(2^128).
Earlier: J8b cache blocking (measured, not adopted), J8b packing by bit-matrix transposes, J2b,
J8e and J8d, J8c, J8, J2a, nimber GF(2) … GF(2^128), ADR-0004, `AGENTS.md`, CI, the move to
MathStruct. Details in git and the Backlog's *Done*.

## Next

Within the approved J8b, only ARM `PMULL` and GPU kernels remain: they need an ARM machine and a
GPU and are not doable in this container. Everything else needs the maintainer's go.

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

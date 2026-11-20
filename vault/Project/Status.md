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
list), and **J8b** (2026-11-20: "Please implement the next item!" again, taken the same way:
J8b was item 1). J8b is a whole section: packing is done, cache blocking measured and not
adopted (below, the third "Please implement the next item!" of 2026-11-20), the rest is open.
Ask before starting a backlog section the maintainer has not approved.

## Done last

- **J8b, cache blocking of the matrix products — measured, not adopted.** First
  `bench/matmul_sizes.jl` (new): the default product from 64 to 2048 does not slow down on this
  VM (GF(2^16) 0.45 → 0.36, GF(2^64) 0.73 → 0.78, GF(2^128) 3.7 → 3.9 ns per multiply-add); its
  260 MiB L3 holds even 64 MiB matrices. Then a blocked `CLMulMatMul` with register tiles
  (1×1 … 4×2): the best tile (2×1) only tied the plain loop at GF(2^16) … GF(2^64) and was
  slower at GF(2^128). The machine code shows the loop is bound by `PCLMULQDQ` throughput
  (about one per term) and the tile's shuffles compete for the same port. So the blocked code
  was **reverted** (it is in git: `79a7ab8`, `20b8cba`); kept are the benchmark, a larger-size
  test (300 × 1100 × 70) and the write-up in `src/bulk.md`, "Larger sizes and cache blocking",
  with three Julia pitfalls met on the way (boxed closure variable, `::Type{T}` unspecialised
  under `@stable`, `@inbounds` not reaching closures).
- **Proposed instead (J8f, not approved):** attack the `PCLMULQDQ` bound: two terms per 64-bit
  product at GF(2^16), `VPCLMULQDQ` (present on this CPU) at every level.
- **Verified:** `Pkg.test()` passes, 732,632 checks, Julia 1.13.1, x86_64-linux, Aqua and JET included. `bench/matmul_sizes.jl` once (the table in `src/bulk.md`); the tile
  comparison one run per tile (noisy, ±10–40%).
- **Not verified:** a CPU with a small L3 (where blocking is expected to matter, not measured);
  Julia 1.10; ARM; CI (not pushed to `full-ci`); the two-terms trick is checked on paper only.

Earlier: J8b packing by bit-matrix transposes (`src/bulk.md`), J2b (the default for `inv`; the
GF(2^32) choice of throughput over latency was confirmed by the maintainer on 2026-11-20), J8e
and J8d, J8c, J8, J2a, nimber GF(2) … GF(2^128), ADR-0004, `AGENTS.md`, CI, the move to
MathStruct. Details in git and the Backlog's *Done*.

## Next

Within the approved J8b (remaining sub-items, in a proposed order):

1. **Bitsliced matrix product** and bitsliced vectors above GF(2^32) (a GF(2^64) prototype was
   slow: tuples of 64 vectors spill).
2. Why inlining both branches of `LazyMulBy`'s `copyto!` slowed the loop (J8d, not inspected).
3. ARM `PMULL` and GPU kernels: need an ARM machine and a GPU; not doable in this container.

Not approved (need the maintainer's go): **J8f** (fewer carry-less products in `CLMulMatMul`;
on the measurements above it is the likelier speed-up for matrix products than item 1);
**FF-P1** (`Problems.md`, `CLMul` at GF(2^64) element by element); **J1a/J1b** (Conway
representation, arbitrary $n$); **J3** (Itoh–Tsujii, extended Euclid as `inv` strategies);
linear algebra on the matrix products (V5).

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

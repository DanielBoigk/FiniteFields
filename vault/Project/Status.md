#overview

> Handoff for the next session: where the work stands, what is next, and what is open. A session
> that starts reads it after `AGENTS.md`; a session that ends updates it. Replace stale text
> rather than appending; history belongs in git.

**Last updated:** 2026-11-19, on branch **`dev-danielb`** (the working branch; pushed).
`plan/agents-handoff` is merged into it and no longer used.

## Stage

Implementation of approved items. Approved and done so far: J0, the nimber fields GF(2) …
GF(2^128), J2a, and J8 (maintainer, 2026-11-19: "please implement J8!"). Ask before starting a
backlog section the maintainer has not approved.

## Done last

- **J8 — bulk and array operations** (`src/bulk.jl`, note `src/bulk.md`), tests and benchmarks
  first:
  - `rand!`/`rand(F, dims)` fill arrays in one call on the storage type: 5–55 times faster than
    per element. Uniform, but different values than before for the same seed.
  - Broadcasting needed no code; tested to fuse without allocations.
  - `MulBy(c)`: multiplication by a constant as byte tables (it is $\mathbb{F}_2$-linear);
    2–11 times faster than `c .* a`, set-up worth 40–240 products.
  - Matrix products: `matmul!(strategy, C, A, B)` with `ScalarMatMul(s)`, `TableMatMul()`,
    `CLMulMatMul()`; `*` and `LinearAlgebra.mul!(C, A, B, α, β)` for strided matrices and vectors.
    `CLMulMatMul` (polynomial basis, unreduced products summed, one reduction per entry) is the
    default from GF(2^16): per multiply-add 0.39 / 0.44 / 0.79 / 3.5 ns at GF(2^16) … GF(2^128),
    against 2.5 / 8.1 / 23 / 79 ns for scalar loops. GF(256) uses `TableMatMul` (0.76 ns).
  - `BitslicedVector` (GF(2) … GF(2^32)): the nimber recursion as an AND/XOR circuit on planes,
    four 64-bit words per SIMD vector; products 4–8 times faster than broadcasting, but packing
    costs 2–6 ns per element.
  - New dependency: `LinearAlgebra` (stdlib).
- **Verified:** `Pkg.test()` passes, 730,160 checks, Julia 1.13.1, x86_64-linux with `PCLMULQDQ`,
  Aqua and JET included. `CLMulMatMul` is also tested with the instruction switched off.
- **Not verified:** Julia 1.10 (`min` in CI); ARM; a CPU without `PCLMULQDQ` (only simulated by the
  flag); any GPU (the claim "Recursive should work in GPU kernels" is an expectation); CI as a
  whole (not pushed to `full-ci`). Benchmarks are a single run on a noisy VM; the decisions rest
  on differences of 6× and more.

Earlier: J2a (`LogTable16`, `CLMul`, default `CLMul` at GF(2^128)), nimber GF(2) … GF(2^128),
ADR-0004, `AGENTS.md`, CI, the move to MathStruct. Details in git.

## Next (proposals, need the maintainer's go)

1. **Defaults for `inv`:** `LogTable16` inverts about 2.3 times faster at GF(2^16) and 1.9 times at
   GF(2^32), but its products have about twice the latency of `MulTable`. A separate default
   for `inv` (`default_inv_strategy`) would take both. Decide whether that API is wanted.
2. **`c .* a` through `MulBy` automatically** for long arrays (2–11 times faster), at the price of
   breaking broadcast fusion for that pattern. Decide whether wanted.
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

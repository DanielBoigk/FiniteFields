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
matrix size!") and **FF-P2** (2026-11-20: "Go ahead with it!"). Ask before starting a backlog
section the maintainer has not approved.

## Done last

- **FF-P2, `CLMulMatMul` for matrix × vector** (maintainer, 2026-11-20: "Go ahead with it!",
  finished after "Yes do so!"). A column-wise path that reads $A$ in storage order without
  transposing it (`clmul_columns_add!`), chosen by `clmul_by_columns(F, m, p, q)`: a vector when
  $A \ge$ 512 KiB, two columns when $A \ge$ 4 MiB and $m, p \ge 256$. The first plan, a rule on
  the number of columns per field (`clmul_column_max`), was dropped: `bench/matvec.jl` showed the
  crossover depends on the size of $A$, and for small $A$ the transposed path is as fast or faster
  even for a vector. $1024 \times 1024$ times a vector: 9.0 / 17.5 / 70 → 1.05 / 2.94 / 14.9 ns
  per multiply-add at GF(2^32) / GF(2^64) / GF(2^128). Two bugs in the work-in-progress commit,
  fixed: `poly_acc_reduce(acc, ::Type{T})` under `@stable` dispatched at runtime per row
  (938 kB allocated per 512 × 512 product), and `poly_acc_zero` was a union for JET. Details in
  `src/bulk.md`, "Matrix × vector".
- **Verified:** `Pkg.test()` passes (733,769 checks, Julia 1.13.1, x86_64-linux);
  `bench/matvec.jl` and `bench/matvec.jl sizes` (single runs); `mul!(y, A, x)` end to end once.
- **Not verified:** CPUs with other cache sizes (the thresholds are this VM's); a CPU without
  `PCLMULQDQ`; other Julia versions; CI. One run at GF(2^64) $128 \times 128$ (5.65 ns, transposed
  path) was slower than the sweep's 3.0–3.1 and was not repeated.

Before that: J8b, why inlining both branches of `LazyMulBy`'s `copyto!` slowed the loop (LLVM's
loop unswitching; `bench/constmul_inline.jl`); J8g (`BitslicedMatMul` the default by size); the
bitsliced matrix product and vectors at GF(2^64)/GF(2^128).
Earlier: J8b cache blocking (measured, not adopted), J8b packing by bit-matrix transposes, J2b,
J8e and J8d, J8c, J8, J2a, nimber GF(2) … GF(2^128), ADR-0004, `AGENTS.md`, CI, the move to
MathStruct. Details in git and the Backlog's *Done*.

## Next

Within the approved J8b, only ARM `PMULL` and GPU kernels remain: they need an ARM machine and a
GPU and are not doable in this container. Everything else needs the maintainer's go.

Not approved (need the maintainer's go): **J8f** (fewer carry-less products in `CLMulMatMul`);
a **cache-blocked transpose** in `CLMulMatMul`'s transposed path (might make the column path
unnecessary and help small $q \ge 3$; an expectation, `src/bulk.md`); faster small bitsliced
products (preallocated temporaries, a `BitslicedMatrix` type keeping the planes; would lower the
thresholds); **FF-P1** (`CLMul` at GF(2^64) element by element); **J1a/J1b** (Conway
representation, arbitrary $n$); **J3** (Itoh–Tsujii, extended Euclid as `inv` strategies);
linear algebra on the matrix products (V5). Seen earlier, not measured: in the scalar GF(2^16)
`LazyMulBy` loop the size of the global `NIM256_MUL` table is reloaded every element.

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

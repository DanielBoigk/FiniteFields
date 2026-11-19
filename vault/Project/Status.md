#overview

> Handoff for the next session: where the work stands, what is next, and what is open. A session
> that starts reads it after `AGENTS.md`; a session that ends updates it. Replace stale text
> rather than appending; history belongs in git.

**Last updated:** 2026-11-19 (second session of the day), on branch `plan/agents-handoff` (made from `main` at `a95d81f`, pushed).

## Stage

**Planning and ideation.** The maintainer has **not** given the go for implementation. No Julia
code exists; Lean has only the one-line `ZMod` smoke test in `lean/FiniteFields.lean`. Work allowed
until the go: notes, plans, specifications, templates, tooling the maintainer asks for.

## Done in the last session (2026-11-19)

- Added `AGENTS.md` (permanent rules and the maintainer's priorities), `CLAUDE.md` (points to it),
  and this handoff.
- **Removed the Julia starter** (`julia/FiniteFields/`, an `Fp{P}` type): it was implemented
  without the maintainer's go. `vault/Specifications/Prime field arithmetic.md` is back to
  *planned*.
- **Fixed `flake.nix`.** The old hooks interpolated `${pkgs.julia-bin}`/`${pkgs.elan}`, which made
  every shell depend on (and download) both. Now:
  - `nix develop` (default) and `.#julia`: use a system Julia if there is one, else build
    `julia-bin` from the pinned nixpkgs on entry. elan is never fetched.
  - `.#lean`: the same for elan (`elan` or `lake` on `PATH` counts as installed). No Julia.
  - `.#quartz`: Git and Node 22 only. `x86_64-darwin` added.
  - Verified here (x86_64-linux, nixpkgs input overridden with a local git checkout of the same
    revision because the session's proxy refused GitHub tarballs): no shell's derivation depends
    on julia or elan; `.#julia` without a system Julia fetched Julia 1.13.1 and had no elan;
    `.#lean` fetched elan 4.2.4 and had no Julia; with a `julia` already on `PATH` the hook used
    it and fetched nothing. **Not verified:** macOS, aarch64, and the normal (non-overridden)
    GitHub fetch of nixpkgs.
  - `flake.lock` is committed. It was written by hand from the content hash of a local export
    of the pinned revision (the proxy blocks GitHub tarballs), and Nix accepted it: it resolved the
    input from the store by that hash and evaluated the shells without the override.
- **`lean/lake-manifest.json` is committed** (no longer git-ignored), written by `lake update`
  with Lean 4.36.0; it pins mathlib to `fc4f59e` (tag `v4.36.0`) and its dependencies.
- **CI:** `lean.yml` (runs `lake build` via `leanprover/lean-action@v1`) and `julia.yml` (tests
  every `julia/<Pkg>/`; skipped while there is none) run on pushes to `main` and on pull requests
  **only when their own paths change**. `full-ci.yml` runs both plus a Quartz build (no deploy) on
  a push to the `full-ci` branch or by hand. `pages.yml` is unchanged. The workflows parse as YAML
  and the action tags exist; **they have not run on GitHub yet.**
- **Cloud environment:** this session moved to a personal environment "FiniteFields" that also
  allows `release.lean-lang.org` and `releases.lean-lang.org` (the Lean toolchain downloads), which
  the Default environment blocked.

## Next

The maintainer wants the arithmetic first (2026-11-19) but has **not yet given the go** to write
code. When it comes, start with backlog J0 (package skeleton, tests, benchmark suite) and J1/J2
for $\mathbb{F}_{2^8}$, in that order, test first. Before that, open questions 1 and 3 below
shape the type design.

## Answers to the maintainer's questions (2026-11-19)

**Can mathlib documentation snippets go into our notes?** Yes, as generated plain Markdown.
Proposal (backlog T1): a small tool that takes a declaration name, asks Lean (against our pinned
mathlib) for its signature, docstring and module, and writes a `> [!mathlib]` callout between
markers (`%% mathlib: Name %%` … `%% /mathlib %%`) into the note, with a link to the mathlib docs
page. Being plain Markdown, it renders in both Obsidian and Quartz; a re-run refreshes it; a
renamed declaration fails loudly. The same tool can sync our own Lean and Julia snippets from
anchors in the source files (like mdBook's `ANCHOR`), so notes never hold stale code. Caveats:
mathlib is Apache-2.0, so quoted docstrings need attribution (a line in the callout); the hosted
mathlib docs track mathlib `master`, not our pin, so a link can point at a declaration that has
moved. A Quartz plugin that fetches docs at build time is possible (a PR to MathStruct/quartz) but
Obsidian would not render it, so the generator is the better first step.

**Can we have mathlib-style API docs for our Lean library, and does that fit with Quartz?** Yes.
mathlib's docs are produced by [doc-gen4](https://github.com/leanprover/doc-gen4), which can run
on any Lake project. It produces a static HTML site that we would publish under a subpath (e.g.
`<site>/lean/`) next to the Quartz vault at the root, the same arrangement Lenticulum.jl uses for
its Documenter API docs under `api/`. MathStruct Quartz already has the fix that lets vault pages
link to same-domain pages outside the vault. Julia API docs (Documenter) would go under
`<site>/julia/` the same way. **Open risk, not checked:** doc-gen4 by default documents every
imported module, i.e. much of mathlib, which is slow and large; whether it can be limited to our
modules (linking out to the hosted mathlib docs) must be tried before we commit to it (T2).

**mathlib facts checked at tag `v4.36.0`** (local sparse checkout, 2026-11-19):
`Mathlib/Algebra/Field/ZMod.lean` (the `Field (ZMod p)` instance);
`Mathlib/FieldTheory/Finite/GaloisField.lean` (`GaloisField p n` as the splitting field of
$X^{p^n}-X$, `GaloisField.card`, `GaloisField.finrank`, `FiniteField.algEquivGaloisField`);
`Mathlib/FieldTheory/Finite/Basic.lean` (`FiniteField.frobeniusAlgHom`, `frobeniusAlgEquiv`);
`Mathlib/FieldTheory/Finite/Trace.lean`; `Mathlib/Algebra/CharP/Two.lean` (`CharTwo.two_eq_zero`,
`CharTwo.neg_one_eq_one`, ...); `Mathlib/Algebra/Polynomial/HasseDeriv.lean`;
`Mathlib/CategoryTheory/Galois/` (Galois categories, fundamental groups). `InnerProductSpace`
(`Mathlib/Analysis/InnerProductSpace/Defs.lean`) requires `[RCLike 𝕜]`, so it cannot be used over
a finite field: the "scalar products over $\mathbb{F}_q$" notes must use bilinear and quadratic
forms (`LinearMap.BilinForm`, `QuadraticForm`) instead. No category of fields (`FieldCat`) was found
under `Mathlib/Algebra/Category/`; if needed it is a full subcategory of `CommRingCat`.

## Decisions by the maintainer (2026-11-19)

- **Milestones:** deferred; the arithmetic is the priority.
- **Field sizes:** first $\mathbb{F}_2$ and the byte-sized fields $\mathbb{F}_{2^8}$,
  $\mathbb{F}_{2^{16}}$, $\mathbb{F}_{2^{32}}$, $\mathbb{F}_{2^{64}}$, $\mathbb{F}_{2^{128}}$.
  Also wanted: a function or macro that constructs $\mathbb{F}_{2^n}$ for **arbitrary $n$**
  (backlog J1b). Asked whether $n = \infty$ is possible: see the answer below (J7).
- **Defining polynomials: both** Conway and standard low-weight ones, one of them the default
  (which one is open question 2).
- **Strategies:** offer both table-based and computed arithmetic where possible; small fields
  probably just lookup tables. Constant-time (cryptographic) use is still under consideration
  (open question 1).
- **Hex** both for `print`/`string` and in the REPL (`show`).
- **Layout:** `julia/FiniteFieldArithmetic.jl/`.
- **`lake-manifest.json`:** commit it (done).
- **CI:** yes, path-filtered, with a branch that runs the full suite (done: `full-ci`).
- **mathlib upstreaming:** decide case by case when it comes up.
- **Site address:** leaning towards `mathstruct.org/FiniteFields` over
  `danielboigk.github.io/FiniteFields`; not final (open question 4).

## Answers to the follow-up questions (2026-11-19)

**Is $n = \infty$ possible?** Yes in a precise sense: the algebraic closure
$\overline{\mathbb{F}}_2 = \bigcup_n \mathbb{F}_{2^n}$. No element needs infinitely many bits;
each lives in some finite $\mathbb{F}_{2^n}$. An implementation stores an element with the degree
of the field it lives in, and an operation on elements of $\mathbb{F}_{2^m}$ and
$\mathbb{F}_{2^n}$ embeds both into $\mathbb{F}_{2^{\operatorname{lcm}(m,n)}}$. This needs a
*compatible* system of embeddings, which is what Conway polynomials provide (or
pseudo-Conway polynomials computed on demand where the tables end). Sizes grow quickly, so it is
an exact-arithmetic tool, not a fast one. In Lean it is mathlib's `AlgebraicClosure (ZMod 2)`. It
is in scope as a later item (J7, V3), after the fixed-size fields.

**What can go wrong when upstreaming to mathlib?** Review takes weeks to months and asks for
mathlib's style, naming and generality, so the result is often restated or moved. A more general
version can clash with existing `RCLike` APIs and need a refactor of mathlib itself. Anything
with `sorry` cannot go. While a PR is open, our copy and mathlib's can diverge, and after a merge
we only get it by bumping our mathlib pin. None of this is harmful; it costs time, which is why
case by case is right.

## Open questions for the maintainer

1. **Constant-time arithmetic** (cryptographic use): in scope or not? If yes, table lookups on
   secret data are out for the affected types, and that splits the strategy design in J2.
2. **Which polynomial is the default:** standard low-weight ones (AES $x^8+x^4+x^3+x+1$, GHASH
   $x^{128}+x^7+x^2+x+1$; interoperable with existing software) or Conway (canonical, compatible
   embeddings, needed for $\overline{\mathbb{F}}_2$)? Recommendation: standard as default for the
   byte-sized fields, Conway selectable and default for the generic constructor.
3. **REPL format:** plain `0x1b`, or `0x1b` plus the field (e.g. `GF(2^8)(0x1b)` in `show`, `0x1b`
   in `print`)? Both are hex; the second shows the type.
4. **Site address:** `mathstruct.org/FiniteFields` needs the repository in the MathStruct
   organisation (or a deploy into it, as Lenticulum.jl does), otherwise GitHub serves it at
   `danielboigk.github.io/FiniteFields`. `site/quartz.config.ts` still has the placeholder.
5. **Milestones** (deferred).

## Known issues

- Wikilinks use two styles (`[[Decisions/…]]` and `[[Project/Decisions/…]]`); pick by-name links
  when the vault is restructured (V0).

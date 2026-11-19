#overview

> Handoff for the next session: where the work stands, what is next, and what is open. A session
> that starts reads it after `AGENTS.md`; a session that ends updates it. Replace stale text
> rather than appending; history belongs in git.

**Last updated:** 2026-11-19, on branch `plan/agents-handoff` (made from `main` at `a95d81f`).

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
  - `flake.lock` is still not committed; the first `nix develop` writes it. Commit it.

## Next

Wait for the maintainer's answers to the open questions below, then turn the plan into backlog
rows that are ready to work on. The plan lives in [[Backlog]] (sections V, T, L, J).

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

## Open questions for the maintainer

1. **Your message ended with an unfinished bullet** ("Please \") in the Quartz section. What did
   you want to add?
2. **Julia layout:** `julia/FiniteFieldArithmetic.jl/` (keeps the current `julia/` folder) or
   Lenticulum's `lib/FiniteFieldArithmetic.jl/`? Recommendation: `julia/`, since this repository
   also has `lean/`.
3. **Field sizes for characteristic 2 first:** which $n$ matter most? Proposal: $\mathbb{F}_{2^8}$
   (AES, Reed–Solomon), $\mathbb{F}_{2^{16}}$, $\mathbb{F}_{2^{64}}$, $\mathbb{F}_{2^{128}}$
   (GHASH), then generic $n \le 128$; anything above 128 bits needs multi-word elements.
4. **Defining polynomials:** Conway polynomials (canonical, compatible embeddings), standard
   low-weight ones (AES's $x^8+x^4+x^3+x+1$), or both, chosen by a type parameter?
5. **Hex printing:** should only `print`/`string` give hex (e.g. `0x1b`), with `show` in the REPL
   giving more context such as `0x1b ∈ GF(2^8)`, or should both be plain hex?
6. **Is cryptographic use in scope** (constant-time arithmetic)? It rules out table lookups for
   secret data and should be decided before the multiplication strategies are designed.
7. **Lean API docs:** try doc-gen4 under a subpath (T2), or are links to the hosted mathlib docs
   plus Lean blocks in the notes enough for now?
8. **mathlib upstreaming:** when we generalise a result mathlib states only for
   $\mathbb{R}$/$\mathbb{C}$, should it also go upstream as a mathlib PR, or stay local?
9. **`lake-manifest.json`** is git-ignored. Recommendation: commit it, so `lake exe cache get`
   and the mathlib revision are reproducible.
10. **CI:** add jobs for `lake build` (and later the Julia tests) beside the Pages workflow?
11. The site's `baseUrl` is the placeholder `example.org/finitefields`. Where will it be hosted
    (e.g. `mathstruct.org/FiniteFields`)?

## Known issues

- Wikilinks use two styles (`[[Decisions/…]]` and `[[Project/Decisions/…]]`); pick by-name links
  when the vault is restructured (V0).

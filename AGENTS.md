# AGENTS.md

Instructions for AI coding agents working in this repository. Humans: everything here applies to
you too. Conventions are adapted from [Lenticulum.jl](https://github.com/MathStruct/Lenticulum.jl)
(its `AGENTS.md` and `vault/Start Here.md`); where this file is silent, follow that repository.

## Read this first

1. [`vault/Project/Status.md`](vault/Project/Status.md) is the handoff: where the work stands,
   what is next, what is open, and which decisions the maintainer has not made yet. Read it after
   this file, and **update it before your session ends** (replace stale text; history belongs in
   git).
2. [`vault/Project/Backlog.md`](vault/Project/Backlog.md) holds the planned work and the ledger
   of open `sorry`s.

**The project is in the planning phase (since 2026-11-19).** Do not write library code (Julia or
Lean beyond the existing smoke test) until the maintainer says implementation may start, and
`Status.md` records it. Until then, work is limited to notes, plans, specifications, templates and
tooling the maintainer asked for. An earlier agent implemented a Julia starter without that go;
it was removed. Do not repeat that.

**Claims need evidence.** A statement in a note is backed by a mathlib declaration, a Lean proof,
a test, a benchmark, or a citation with an exact theorem number. Anything else is labelled a
conjecture or an expectation. When you find something wrong, fix it or record it in
`vault/Project/Problems.md`; never work around it silently.

## The repository

| where | what |
|---|---|
| `vault/` | the Obsidian vault, published with Quartz: concept notes, specifications, project records |
| `vault/Project/` | `Status.md` (handoff), `Backlog.md`, `Problems.md`, `Decisions/` (ADRs) |
| `vault/templates/` | note templates; start every new note from one |
| `lean/` | Lean 4 + mathlib (`lean-toolchain` and the mathlib tag in `lakefile.lean` move together) |
| `julia/` | Julia packages, one directory each (`julia/FiniteFieldArithmetic.jl/` first); does not exist yet |
| `site/` | the Quartz build: `build.sh`, `quartz.config.ts`, `quartz.layout.ts` |
| `flake.nix`, `flake.lock` | dev shells; see the README |
| `.github/workflows/` | CI, see below |

## Priorities set by the maintainer

- **Purpose:** the general theory of finite fields, cryptanalysis perhaps, and machine-learning
  "learners" built on finite fields (the maintainer's experiment). Implementing cryptographic
  algorithms is not a goal, and neither is constant-time arithmetic.

- **Characteristic 2 first.** Fields $\mathbb{F}_{2^n}$ come before odd characteristic, in the
  vault, in Lean and in Julia.
- **Categorical view where possible.** Notes state a concept categorically when that is natural,
  or say how it relates to category theory; Lean adds category-theoretic definitions and theorems
  (mathlib's `CategoryTheory`) beside the classical ones.
- **Contrast with $\mathbb{R}$ and $\mathbb{C}$.** Where finite-field theory differs (inner
  products, orderings, absolute values, differentiability, Fourier analysis, ...), there is a note
  saying exactly what survives and what fails.

## Vault conventions

- **One concept per note.** Small, granular notes, each explaining exactly one definition,
  theorem, proof, example, algorithm or comparison. Split a note that grows a second subject.
- **Tag line first**, from this vocabulary (several allowed): `#definition`, `#theorem`,
  `#proof`, `#example`, `#algorithm`, `#construction`, `#comparison`, `#categorical`,
  `#characteristic-2`, `#implementation`, `#design`, `#open-problem`, `#overview`, `#reference`.
  Extend the vocabulary here, not ad hoc.
- **Then a summary** blockquote (two or three lines), then a **sources** block:
  `> mathlib:` declarations with links, `> Lean:` our files/declarations, `> Julia:` our code and
  API docs, `> Literature:` with exact theorem numbers.
- **Lean formalization in the note.** A note whose subject is formalized shows the Lean statement
  (mathlib's or ours) in a ` ```lean ` block, and links the mathlib docs page
  (`https://leanprover-community.github.io/mathlib4_docs/<Module/Path>.html#<Decl.Name>`).
  Lean and Julia versions of the same thing may sit side by side in a ` ```tabs ` block
  (MathStruct Quartz plugin).
- Math in LaTeX; display `$$` fences on their own lines. Diagrams in ` ```tikz ` blocks.
  Wikilinks by note name.
- Quartz: use [MathStruct/quartz](https://github.com/MathStruct/quartz), pinned in
  `site/build.sh`. A missing Quartz feature is a pull request there (ask the maintainer first),
  not a local patch.

## Lean conventions

- **mathlib first.** Search mathlib before defining anything; import the narrowest module; cite
  the exact declaration and link it. Policy details: `vault/Project/Lean and mathlib.md`.
- **Generalise what mathlib over-specialises.** A result stated in mathlib only for
  $\mathbb{R}$, $\mathbb{C}$ or `RCLike` that holds over any field (or any field of a given
  characteristic) is restated and proved generically in our library, with a link to the
  specialised original.
- **Statements first, proofs second.** A definition or theorem may land with `sorry` in a first
  sweep. **Every `sorry` gets a row in the Backlog's sorry ledger** (file, declaration, note) in
  the same change; the proofs come in a later sweep. Never leave an unlisted `sorry`.

## Julia conventions (for when implementation starts)

- Orient on [Lenticulum.jl](https://github.com/MathStruct/Lenticulum.jl): one environment per
  package, `test/aqua.jl` closing every suite, a `bench/` environment, implementation notes
  (`x.md` beside `src/x.jl`).
- **Test first.** Write the test (and, for a performance-relevant function, the benchmark) before
  the implementation. Profile performance-critical code (`Profile`, `@profview`, allocation
  profiler) and record the numbers in a note together with the script that produced them.
- **Quality gates:** Aqua.jl, JET.jl (`@test_opt`/`report_package`), DispatchDoctor.jl
  (`@stable` on hot paths), and similar tools in every test suite.
- **Multiple dispatch over needless specialisation.** Write generic methods; specialise only where
  a benchmark shows it pays.
- **No integers wider than necessary.** Element storage uses the smallest unsigned type that
  fits (`UInt8` for $\mathbb{F}_{2^8}$, ...). Offer several arithmetic strategies where they
  trade off (tables, carry-less multiplication, ...), selected by dispatch.
- **Characteristic-2 elements print in hexadecimal**, both with `print` and in the REPL (`show`).
- The first package is **`FiniteFieldArithmetic.jl`**: the basic field operations.

## CI

`lean.yml` and `julia.yml` run on pushes to `main` and on pull requests, each **only when its own
paths change**, so a vault or Julia change does not rebuild against mathlib. To run everything
(Lean, every Julia package, a Quartz build) push to the `full-ci` branch
(`git push origin HEAD:full-ci --force`) or start "Full CI" from the Actions tab; do that before
merging a change that spans components. `pages.yml` deploys the vault on pushes to `main`.

## Commits and handoff

A task ends with `Status.md` and the Backlog updated. Commit messages say what changed, why, and
what was **not** checked. Do not push or open pull requests unless asked.

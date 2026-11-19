# FiniteFields

> **Planning phase.** No library code is written yet; see [`AGENTS.md`](AGENTS.md) and the
> handoff in [`vault/Project/Status.md`](vault/Project/Status.md).

A small, proof-oriented monorepo for finite-field mathematics. It brings together:

- an Obsidian vault for documentation, decisions, problems, backlog, and specifications;
- Julia packages, beginning with `FiniteFieldArithmetic.jl` (fields of characteristic 2 first)
  and growing toward linear algebra and analysis;
- Lean 4 formalizations, reusing mathlib whenever it already contains the definition or proof.

The vault is published with [MathStruct's Quartz](https://github.com/MathStruct/quartz), pinned
to a commit that includes its TikZ, Typst/LaTeX, and Obsidian plugins. Project notes follow the
decision/problem/backlog separation used by
[Lenticulum.jl](https://github.com/MathStruct/Lenticulum.jl).

## Development environment

Install [Nix with flakes enabled](https://nixos.org/download/), then pick the shell for the
language you are working on:

```sh
nix develop           # default: Julia plus Node.js 22 for the site; no Lean tooling
nix develop .#julia   # Julia only
nix develop .#lean    # Lean only (elan); no Julia
nix develop .#quartz  # vault/site work: Node.js 22 and Git
```

Toolchains are fetched lazily. On entering a shell, the hook looks for the tool on your system
(`julia`; `elan` or `lake`) and uses it if present. Only when it is missing does the hook build
the Nix copy from the pinned nixpkgs and put it on `PATH`, so a system install is never shadowed
or downloaded twice, and elan is only ever fetched by `.#lean`. `elan` (or an existing Lean)
installs the Lean version pinned in [`lean/lean-toolchain`](lean/lean-toolchain).

## Julia

Nothing is implemented yet. The first package will be `FiniteFieldArithmetic.jl`: the basic field
operations, with fields of characteristic 2 first. The plan is in the
[backlog](vault/Project/Backlog.md) (section J) and the conventions in [`AGENTS.md`](AGENTS.md).

## Lean 4

```sh
cd lean
lake exe cache get
lake build
```

[`FiniteFields.lean`](lean/FiniteFields.lean) imports mathlib's `ZMod` field instance instead of
redefining it. See [Lean and mathlib](vault/Project/Lean%20and%20mathlib.md) for the reuse policy.

## Documentation vault

Open [`vault/`](vault/) as an Obsidian vault. Start at
[`Start Here.md`](vault/Start%20Here.md). To build or preview the Quartz site:

```sh
site/build.sh
site/build.sh --serve
```

The first build checks out the pinned Quartz revision into `site/.quartz/` and installs its
locked npm dependencies. Generated files and the Quartz checkout are ignored by Git.

On pushes to `main` that touch `vault/` or `site/`, the
[publish workflow](.github/workflows/pages.yml) rebuilds the site and deploys it to the
`gh-pages` branch.

## Project workflow

1. Write user-visible behavior and acceptance criteria in `vault/Specifications/`.
2. Record consequential design choices in `vault/Project/Decisions/`.
3. Track discovered defects in `vault/Project/Problems.md` and planned work in
   `vault/Project/Backlog.md`.
4. Link Julia implementations and Lean formalizations to their specifications. Prefer a mathlib
   definition or theorem when it already covers the result; document the exact import/name.

This is a starter repository, not a claim that every finite-field construction or theorem has
already been implemented or formalized.

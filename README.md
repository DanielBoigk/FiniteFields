# FiniteFields

A small, proof-oriented monorepo for finite-field mathematics. It brings together:

- an Obsidian vault for documentation, decisions, problems, backlog, and specifications;
- Julia packages, beginning with arithmetic in prime fields and growing toward linear algebra
  and analysis;
- Lean 4 formalizations, reusing mathlib whenever it already contains the definition or proof.

The vault is published with [MathStruct's Quartz](https://github.com/MathStruct/quartz), pinned
to a commit that includes its TikZ, Typst/LaTeX, and Obsidian plugins. Project notes follow the
decision/problem/backlog separation used by
[Lenticulum.jl](https://github.com/MathStruct/Lenticulum.jl).

## Development environment

Install [Nix with flakes enabled](https://nixos.org/download/), then run:

```sh
nix develop
```

The shell provides Julia, Node.js 22, Git, and `elan`. `elan` installs the Lean version pinned in
[`lean/lean-toolchain`](lean/lean-toolchain).

## Julia

```sh
julia --project=julia/FiniteFields -e 'using Pkg; Pkg.test()'
```

The starter package implements `Fp{p}` for prime `p`, with checked construction and basic field
arithmetic. It deliberately starts with prime fields; polynomial-basis extension fields are
tracked in the vault backlog.

## Lean 4

```sh
cd lean
lake exe cache get
lake build
```

[`FiniteFields.lean`](lean/FiniteFields.lean) imports mathlib's `ZMod` field instance instead of
redefining it. See [Lean and mathlib](vault/Project/Lean and mathlib.md) for the reuse policy.

## Documentation vault

Open [`vault/`](vault/) as an Obsidian vault. Start at
[`Start Here.md`](vault/Start%20Here.md). To build or preview the Quartz site:

```sh
site/build.sh
site/build.sh --serve
```

The first build checks out the pinned Quartz revision into `site/.quartz/` and installs its
locked npm dependencies. Generated files and the Quartz checkout are ignored by Git.

## Project workflow

1. Write user-visible behavior and acceptance criteria in `vault/Specifications/`.
2. Record consequential design choices in `vault/Project/Decisions/`.
3. Track discovered defects in `vault/Project/Problems.md` and planned work in
   `vault/Project/Backlog.md`.
4. Link Julia implementations and Lean proofs to their specifications. Prefer a mathlib
   definition or theorem when it already covers the result; document the exact import/name.

This is a starter repository, not a claim that every finite-field construction or theorem has
already been implemented or formalized.

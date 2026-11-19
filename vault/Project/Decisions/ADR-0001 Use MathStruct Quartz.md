# ADR-0001: Use MathStruct Quartz

**Status:** accepted

## Context

The project needs a browsable Obsidian vault with mathematical documentation and TikZ diagrams.

## Decision

Build the vault with [MathStruct/quartz](https://github.com/MathStruct/quartz), pinned to
`6c85c0950ac0bd9567d2509dd1d50d38e4e6696c`. Its MathStruct preset supplies TikZ and Typst/LaTeX
support as well as Obsidian-oriented plugins. The wrapper in `site/build.sh` fetches that exact
commit.

## Consequences

- Site builds are reproducible against a known Quartz revision.
- Updating Quartz is an explicit change to the pin and should be tested against the vault.
- The upstream checkout and generated site are disposable build outputs, not project sources.

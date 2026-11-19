# ADR-0002: Reuse mathlib

**Status:** accepted

## Context

Lean is intended to state and prove the project's mathematics, but finite-field constructions and
theorems may already exist in mathlib.

## Decision

Search mathlib before adding a definition or proof. For prime fields, use `ZMod p` and its `Field`
instance from `Mathlib.Algebra.Field.ZMod`. New project-specific results should be proved using
those existing structures and should link to the relevant specification.

## Consequences

The Lean starter depends on mathlib at the same release as its pinned Lean toolchain. Any new
formalization should say whether it is new or a direct reuse of a mathlib declaration.

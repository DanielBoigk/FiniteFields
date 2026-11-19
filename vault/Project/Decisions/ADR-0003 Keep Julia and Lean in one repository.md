# ADR-0003: Keep Julia and Lean in one repository

**Status:** accepted

## Context

Julia is the executable library and Lean is its formal specification/proof environment. Splitting
them would make it harder to keep mathematical requirements, implementations, and proofs aligned.

## Decision

Keep the Julia package, Lean project, and specification vault together. Link related artifacts
from their specifications and track work that spans components as one backlog item.

## Consequences

Julia dependencies and Lean dependencies remain independently managed. The Nix development shell
provides their toolchains without replacing their native project manifests.

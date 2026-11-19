# FiniteFields

This vault is the project record for a finite-fields library with Julia implementations and Lean
formalizations.

Agents and contributors: read `AGENTS.md` in the repository root, then [[Project/Status]].

## Project areas

- [[Project/Status]] — handoff: where the work stands, what is next, open questions
- [[Project/Decisions]] — decision index and rationale
- [[Project/Problems]] — known gaps and defects
- [[Project/Backlog]] — planned work
- [[Specifications/Prime field arithmetic]] — first behavior specification
- [[Project/Lean and mathlib]] — formalization and reuse policy

## Working rule

Keep the mathematical statement independent of its implementation. Link a specification to the
Julia code and Lean theorem that realize it; when mathlib already supplies a definition or proof,
cite the import and declaration instead of creating a competing copy.

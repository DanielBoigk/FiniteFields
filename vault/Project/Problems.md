# Problems

Record each discovered defect with an identifier, a minimal
reproduction, expected and actual behavior, and links to the affected specification and code.

Do not use this log as a substitute for feature planning; planned work belongs in
[[Backlog|the backlog]].

## FF-P1 — `CLMul` at GF(2^64) loses its speed in the element wrapper (open, 2026-11-19)

- **Where:** `julia/FiniteFieldArithmetic.jl/src/nimber_fast.jl`, `nim_mul_bits(::CLMul, …, Val(6))`
  and `mul(CLMul(), a::GF2_64, b::GF2_64)`; measurements in `src/nimber_fast.md`.
- **Reproduction:** broadcast over 4096 random pairs, in one process. On `UInt64` bit patterns,
  `nim_mul_bits(CLMul(), x, y, Val(6))` takes 10.4–11.7 ns per product against 18–20 ns for
  `MulTable`. Wrapped as `GF2_64` (`mul(CLMul(), a, b)` or
  `unsafe_frombits(GF2_64, nim_mul_bits(…))`) it takes 17–19 ns. The native code of both
  functions is nearly identical (168 vs 179 instructions, the same single out-of-line call).
- **Expected:** the wrapper (an `isbits` struct around the `UInt64`) costs nothing.
- **Actual:** about 7 ns per product in some runs. The VM's run-to-run noise (up to about 40%)
  makes it unclear how much is real: one run of the bare variant also took 18 ns.
- **Consequence:** GF(2^64) keeps `MulTable` as default. Investigate on a quiet machine
  (`code_llvm` of the whole broadcast loop, `perf`), then reconsider the default.

## FF-P2 — `CLMulMatMul` is slow for matrix × vector (open, 2026-11-20)

- **Where:** `julia/FiniteFieldArithmetic.jl/src/bulk.jl`, `matmul_add!(::CLMulMatMul, …)`, the
  default for `A * x` above GF(256) (below the bitsliced threshold, which a vector never reaches).
- **Reproduction:** `bench/matmul_threshold.jl`, rows `… × 1`: $1024 \times 1024$ times a vector
  takes 9.0 ns per multiply-add at GF(2^32), 17.5 at GF(2^64), 70 at GF(2^128), against 0.3, 0.6
  and 3 per multiply-add for square products.
- **Expected:** about the cost of a square product per multiply-add.
- **Actual:** the conversion of all of $A$ to the polynomial basis (and $A$ transposed) is paid
  for one column; for a vector it costs more than the products. (An explanation from the code,
  not profiled.)
- **Consequence:** none for correctness. Options: a matrix × vector path without converting $A$
  (scalar `CLMul` products, or `MulBy`-style tables of the vector entries), chosen by $q$.

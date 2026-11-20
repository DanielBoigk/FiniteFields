"""
    FiniteFieldArithmetic

Arithmetic in finite fields, characteristic 2 first. Implemented so far: the nimber fields
``\\mathbb{F}_2 \\subset \\mathbb{F}_4 \\subset \\dots \\subset \\mathbb{F}_{2^{128}}``
([`GF2`](@ref), [`GF4`](@ref), [`GF16`](@ref), [`GF256`](@ref), [`GF2_16`](@ref),
[`GF2_32`](@ref), [`GF2_64`](@ref), [`GF2_128`](@ref)), nested as integers.

Multiplication strategies ([`MulStrategy`](@ref): [`Recursive`](@ref), [`MulTable`](@ref),
[`LogTable`](@ref), [`LogTable16`](@ref), [`CLMul`](@ref)) are reached through [`mul`](@ref) and
[`inverse`](@ref); `*` uses [`default_strategy`](@ref), `inv` (and `/`)
[`default_inv_strategy`](@ref).

Bulk operations (`src/bulk.md`): fast `rand!` for arrays, [`MulBy`](@ref) (multiplication by a
constant), matrix products ([`matmul!`](@ref) with a [`MatMulStrategy`](@ref), among them
[`BitslicedMatMul`](@ref) (`src/bitsliced_matmul.md`); `*` and
`LinearAlgebra.mul!` use [`default_matmul_strategy`](@ref)), [`BitslicedVector`](@ref), and the
strategies of `c .* a` ([`ConstMulStrategy`](@ref), chosen by [`default_constmul_strategy`](@ref)).

The smallest field of an element (`src/subfield.md`): [`minimal_level`](@ref),
[`minimal_degree`](@ref), [`reduce_field`](@ref).

Type stability of the arithmetic is enforced by DispatchDoctor's `@stable` (mode `"error"`); set
the preference `instability_check = "disable"` for this package to switch the checks off.
"""
module FiniteFieldArithmetic

using DispatchDoctor: @stable
using LinearAlgebra: LinearAlgebra
using Random: Random

export BinaryFieldElement, Nimber, GF2, GF4, GF16, GF256, GF2_16, GF2_32, GF2_64, GF2_128
export characteristic, degree, elements
export minimal_level, minimal_degree, reduce_field
export MulBy, BitslicedVector

include("binaryfield.jl")
include("nimber.jl")
include("nimber_fast.jl")
include("bulk.jl")
include("bitsliced_matmul.jl")
include("subfield.jl")

end

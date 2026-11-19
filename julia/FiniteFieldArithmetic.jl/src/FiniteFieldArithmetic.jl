"""
    FiniteFieldArithmetic

Arithmetic in finite fields, characteristic 2 first. Implemented so far: the nimber fields
``\\mathbb{F}_2 \\subset \\mathbb{F}_4 \\subset \\mathbb{F}_{16} \\subset \\mathbb{F}_{256}``
([`GF2`](@ref), [`GF4`](@ref), [`GF16`](@ref), [`GF256`](@ref)), nested as integers.

Multiplication strategies ([`MulStrategy`](@ref): [`Recursive`](@ref), [`MulTable`](@ref),
[`LogTable`](@ref)) are reached through [`mul`](@ref) and [`inverse`](@ref); `*` and `inv` use
[`default_strategy`](@ref).

Type stability of the arithmetic is enforced by DispatchDoctor's `@stable` (mode `"error"`); set
the preference `instability_check = "disable"` for this package to switch the checks off.
"""
module FiniteFieldArithmetic

using DispatchDoctor: @stable
using Random: Random

export BinaryFieldElement, Nimber, GF2, GF4, GF16, GF256
export characteristic, degree, elements

include("binaryfield.jl")
include("nimber.jl")

end

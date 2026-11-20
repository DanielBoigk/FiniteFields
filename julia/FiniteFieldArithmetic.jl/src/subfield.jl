# J9: the smallest field of an element. Theory and the nimber shortcut: subfield.md beside this
# file. Other representations (Conway, polynomial bases, J1a/J1b; the algebraic closure, J7) add
# methods later, through the Frobenius: the smallest d with x^(2^d) = x.

"""
    minimal_degree(x::BinaryFieldElement) -> Int

The degree ``[\\mathbb{F}_2(x) : \\mathbb{F}_2]`` of the smallest subfield containing `x`, which
is also the degree of the minimal polynomial of `x` over ``\\mathbb{F}_2`` and the smallest
``d \\ge 1`` with ``x^{2^d} = x``. It divides `degree(typeof(x))`.

For nimbers it is ``2^j`` with `j = minimal_level(x)`, read off the bit length; no field
arithmetic is done.

```jldoctest
julia> minimal_degree(GF2_64(0x1b))
8
```
"""
function minimal_degree end

"""
    minimal_level(x::Nimber) -> Int

The smallest level ``j`` with ``x \\in \\mathbb{F}_{2^{2^j}}``, i.e. the smallest ``j`` with
``x < 2^{2^j}`` as an integer (the nimber fields are nested as integers, and every subfield of
``\\mathbb{F}_{2^{2^K}}`` is one of them). `0` for `0` and `1`.
"""
function minimal_level end

"""
    reduce_field(x::Nimber)

`x` as an element of the smallest field containing it: the nimber level
[`minimal_level`](@ref)`(x)`, same value, `== x`. `reduce_field(GF2_128(3)) === GF4(3)`.

**Type-unstable by design** (maintainer, 2026-11-20): the result type depends on the value, so
it is not covered by DispatchDoctor's `@stable`. Use it at the boundary of a computation, not in
a hot loop; [`minimal_level`](@ref) and [`minimal_degree`](@ref) are type-stable.
"""
function reduce_field end

@stable default_mode = "error" begin
    function minimal_level(x::Nimber{K, T}) where {K, T}
        b = bits(x)
        len = 8 * sizeof(T) - leading_zeros(b)                 # bit length, 0 for x = 0
        # the smallest j with len ≤ 2^j: ⌈log₂ len⌉, and 0 for len ≤ 1
        return 8 * sizeof(Int) - leading_zeros(max(len, 1) - 1)
    end

    minimal_degree(x::Nimber) = 1 << minimal_level(x)
end

# One branch per level, each with a concrete type, instead of building `Nimber{j, S}` at
# runtime: the result is a Union of at most K + 1 types, not a dynamically constructed type
# (bench/subfield.jl, random elements: 180–205 ns per element for the latter, 3.5–4.7 ns for
# this).
function reduce_field(x::Nimber{K, T}) where {K, T}
    b = bits(x)
    j = minimal_level(x)
    j <= 0 && return unsafe_frombits(GF2, b % UInt8)
    j == 1 && return unsafe_frombits(GF4, b % UInt8)
    j == 2 && return unsafe_frombits(GF16, b % UInt8)
    j == 3 && return unsafe_frombits(GF256, b % UInt8)
    j == 4 && return unsafe_frombits(GF2_16, b % UInt16)
    j == 5 && return unsafe_frombits(GF2_32, b % UInt32)
    j == 6 && return unsafe_frombits(GF2_64, b % UInt64)
    return unsafe_frombits(GF2_128, b % UInt128)
end

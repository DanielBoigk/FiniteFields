# What every representation of a field of characteristic 2 shares: the element is a bit pattern,
# addition is XOR, negation is the identity, printing is hexadecimal. A representation
# (nimbers now; Conway and standard polynomial bases later) subtypes `BinaryFieldElement` and
# provides `bits`, `storagetype`, `degree`, `unsafe_frombits`, `*` and `inv`.

"""
    BinaryFieldElement <: Number

Supertype of elements of finite fields of characteristic 2, ``\\mathbb{F}_{2^n}``, stored as a
bit pattern in an unsigned integer. Subtypes provide [`bits`](@ref), [`degree`](@ref),
`storagetype`, `unsafe_frombits`, `*` and `inv`; addition (XOR), subtraction, negation,
equality, hashing, powers, division and hexadecimal printing are defined here once for all of
them.
"""
abstract type BinaryFieldElement <: Number end

"""
    bits(x::BinaryFieldElement)

The bit pattern of `x`, an unsigned integer of the field's storage type.
"""
function bits end

"""
    degree(F::Type{<:BinaryFieldElement})

The degree ``n`` of ``F = \\mathbb{F}_{2^n}`` over ``\\mathbb{F}_2``.
"""
function degree end

"""
    characteristic(F::Type{<:BinaryFieldElement}) = 2
"""
characteristic(::Type{<:BinaryFieldElement}) = 2

"""
    elements(F::Type{<:BinaryFieldElement})

All ``2^n`` elements of `F`, in the order of their bit patterns.
"""
elements(::Type{F}) where {F <: BinaryFieldElement} =
    (unsafe_frombits(F, storagetype(F)(i)) for i in 0:(2^degree(F) - 1))

# Whether the bit pattern fills the storage type (GF(256) in a UInt8 does, GF(16) does not).
fills_storage(::Type{F}) where {F <: BinaryFieldElement} = degree(F) == 8 * sizeof(storagetype(F))

# Order of the multiplicative group, as an Int. Valid while degree(F) ≤ 62, which covers every
# field implemented so far; wider fields need a wider exponent type here.
multiplicative_order(::Type{F}) where {F <: BinaryFieldElement} = (1 << degree(F)) - 1

@stable default_mode = "error" begin
    Base.:+(x::F, y::F) where {F <: BinaryFieldElement} = unsafe_frombits(F, bits(x) ⊻ bits(y))
    Base.:-(x::F, y::F) where {F <: BinaryFieldElement} = x + y
    Base.:-(x::BinaryFieldElement) = x
    Base.conj(x::BinaryFieldElement) = x
    Base.zero(::Type{F}) where {F <: BinaryFieldElement} = unsafe_frombits(F, zero(storagetype(F)))
    Base.one(::Type{F}) where {F <: BinaryFieldElement} = unsafe_frombits(F, one(storagetype(F)))
    Base.iszero(x::BinaryFieldElement) = iszero(bits(x))
    Base.isone(x::BinaryFieldElement) = isone(bits(x))
    Base.:(==)(x::F, y::F) where {F <: BinaryFieldElement} = bits(x) == bits(y)
    Base.:/(x::F, y::F) where {F <: BinaryFieldElement} = x * inv(y)

    function Base.:^(x::BinaryFieldElement, n::Integer)
        if iszero(x)
            n < 0 && throw(DivideError())
            return n == 0 ? one(x) : x
        end
        # x^(2^n - 1) = 1 for x ≠ 0, so the exponent only matters modulo the group order;
        # reducing it also turns a negative exponent into a positive one.
        e = Int(mod(n, multiplicative_order(typeof(x))))
        return Base.power_by_squaring(x, e)
    end
end

# Elements of different types are compared in their common promotion, if there is one (nimber
# levels: the larger level, since the subfields are nested). Unrelated representations promote to
# no concrete type and are never equal: there is no implicit isomorphism between them.
function Base.:(==)(x::BinaryFieldElement, y::BinaryFieldElement)
    P = promote_type(typeof(x), typeof(y))
    isconcretetype(P) || return false
    return convert(P, x) == convert(P, y)
end

# A field element is never equal to an integer, a float or any other number: there is no ring
# map that would make the comparison meaningful in general (an integer could mean a bit pattern
# or the image of ℤ). Returning false rather than throwing keeps `in`, `Set` and `Dict` working
# on mixed collections.
Base.:(==)(::BinaryFieldElement, ::Number) = false
Base.:(==)(::Number, ::BinaryFieldElement) = false

# Hash the bit pattern only: representations whose subfields are nested as integers (nimbers)
# compare equal across levels, so their hashes must agree too.
Base.hash(x::BinaryFieldElement, h::UInt) = hash(bits(x), hash(:BinaryFieldElement, h))

# Uniform random elements.
function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{F}) where {F <: BinaryFieldElement}
    S = storagetype(F)
    mask = fills_storage(F) ? typemax(S) : (one(S) << degree(F)) - one(S)
    return unsafe_frombits(F, rand(rng, S) & mask)
end

# Printing: hexadecimal, one digit per four bits of the field (maintainer, 2026-11-19).
function Base.print(io::IO, x::BinaryFieldElement)
    n = degree(typeof(x))
    print(io, "0x", string(bits(x); base = 16, pad = cld(n, 4)))
    return nothing
end

# In the REPL the field is named when the hex value alone does not determine it (the field does
# not fill its storage), unless the context already states the element type (typed arrays).
function Base.show(io::IO, x::BinaryFieldElement)
    F = typeof(x)
    print(io, x)
    if !fills_storage(F) && get(io, :typeinfo, Any) !== F
        n = degree(F)
        print(io, n == 1 ? " ∈ GF(2)" : " ∈ GF(2^$n)")
    end
    return nothing
end

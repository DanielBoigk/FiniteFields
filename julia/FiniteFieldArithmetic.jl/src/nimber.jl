# Nimber fields GF(2^(2^K)). Theory, algorithms and measurements: nimber.md beside this file.

"""
    Nimber{K, T<:Unsigned} <: BinaryFieldElement

An element of the nimber field ``\\mathbb{F}_{2^{2^K}}``: the integers ``0 \\le a < 2^{2^K}``
with nim-addition (XOR) and nim-multiplication (Conway, *On Numbers and Games*, ch. 6), stored in
`T`, the smallest unsigned type that fits (`UInt8` for ``K \\le 3``).

The nimber fields are nested *as integers*: ``\\mathbb{F}_{2^{2^K}}`` consists of the integers
below ``2^{2^K}``, so the embedding into a larger nimber field keeps the bit pattern. Mixed
arithmetic promotes to the larger field accordingly, and elements of different levels with the
same value compare equal.

The default representation of the fields ``\\mathbb{F}_{2^{2^K}}`` (maintainer, 2026-11-19).
Implemented so far: ``K \\le 3``, i.e. [`GF2`](@ref), [`GF4`](@ref), [`GF16`](@ref),
[`GF256`](@ref).

Construct from the integer value (= the bit pattern): `GF256(0x1b)`, `GF256(27)`.
"""
struct Nimber{K, T <: Unsigned} <: BinaryFieldElement
    bits::T
    function Nimber{K, T}(b::T, ::Val{:unchecked}) where {K, T <: Unsigned}
        nimber_storage(K) === T || throw(ArgumentError(
            "Nimber{$K, $T}: GF(2^$(2^K)) is stored in $(nimber_storage(K))"))
        return new{K, T}(b)
    end
end

"""GF(2) as nimbers: `Nimber{0, UInt8}`."""
const GF2 = Nimber{0, UInt8}
"""GF(4) as nimbers: `Nimber{1, UInt8}`; contains [`GF2`](@ref) as the values 0, 1."""
const GF4 = Nimber{1, UInt8}
"""GF(16) as nimbers: `Nimber{2, UInt8}`; contains [`GF4`](@ref) as the values 0–3."""
const GF16 = Nimber{2, UInt8}
"""GF(256) as nimbers: `Nimber{3, UInt8}`; contains [`GF16`](@ref) as the values 0–15."""
const GF256 = Nimber{3, UInt8}

# The storage type of level K. Levels above 3 (GF(2^16) and up) are not implemented yet.
function nimber_storage(K::Int)
    0 <= K <= 3 && return UInt8
    throw(ArgumentError("nimber fields above GF(2^8) (level $K) are not implemented yet"))
end

bits(x::Nimber) = x.bits
storagetype(::Type{Nimber{K, T}}) where {K, T} = T
degree(::Type{Nimber{K, T}}) where {K, T} = 1 << K
unsafe_frombits(::Type{Nimber{K, T}}, b::T) where {K, T} = Nimber{K, T}(b, Val(:unchecked))

function Nimber{K, T}(n::Integer) where {K, T <: Unsigned}
    0 <= n < (1 << (1 << K)) ||
        throw(ArgumentError("$n is not an element of GF(2^$(1 << K)) (needs 0 ≤ n < 2^$(1 << K))"))
    return unsafe_frombits(Nimber{K, T}, T(n))
end

# Embedding into a larger level keeps the bits; restriction to a smaller one checks they fit.
Nimber{K, T}(x::Nimber) where {K, T <: Unsigned} = Nimber{K, T}(bits(x))

Base.promote_rule(::Type{Nimber{K1, T1}}, ::Type{Nimber{K2, T2}}) where {K1, T1, K2, T2} =
    K1 >= K2 ? Nimber{K1, T1} : Nimber{K2, T2}

# --- Strategies -----------------------------------------------------------------------------

"""
    MulStrategy

How multiplication and inversion are computed. Every strategy gives the same results; they
differ in speed and memory, measured in `bench/` and recorded in `src/nimber.md`.

- [`Recursive`](@ref): Karatsuba-style recursion over the levels, no tables.
- [`MulTable`](@ref): a full 256 × 256 product table (64 KiB) and a 256-entry inverse table.
- [`LogTable`](@ref): discrete logarithm and exponential tables for a generator of
  ``\\mathbb{F}_{256}^\\times`` (768 bytes).

The table strategies serve every level up to ``\\mathbb{F}_{256}``: since the subfields are
nested as integers, the tables of ``\\mathbb{F}_{256}`` restricted to the values below
``2^{2^K}`` are those of ``\\mathbb{F}_{2^{2^K}}``.
"""
abstract type MulStrategy end

"""
    Recursive()

Nim-multiplication by recursion over the levels. With ``h = 2^{K-1}``, ``F = 2^h`` and
``a = a_1 F + a_0``, ``b = b_1 F + b_0`` (halves in the subfield of level ``K-1``):
``a \\otimes b = (e \\oplus d)\\,F + (d \\oplus c \\otimes F/2)`` with ``c = a_1 \\otimes b_1``,
``d = a_0 \\otimes b_0``, ``e = (a_0 \\oplus a_1) \\otimes (b_0 \\oplus b_1)``, because
``F \\otimes F = F \\oplus F/2``; the product with the constant ``F/2`` has its own, cheaper
recursion. Inversion uses the norm to the subfield.
No memory beyond registers; the reference implementation for the table strategies.
"""
struct Recursive <: MulStrategy end

"""
    MulTable()

Lookup in a full product table of ``\\mathbb{F}_{256}`` (64 KiB) and an inverse table.
"""
struct MulTable <: MulStrategy end

"""
    LogTable()

``a \\otimes b = \\exp(\\log a + \\log b)`` with logarithm and exponential tables for a fixed
generator of ``\\mathbb{F}_{256}^\\times``; zero handled separately.
"""
struct LogTable <: MulStrategy end

"""
    default_strategy(F)

The strategy behind `*` and `inv` for `F`, chosen by the benchmarks in `bench/`
(`src/nimber.md`, "Measurements").
"""
default_strategy(::Type{<:Nimber}) = MulTable()

# --- Recursive ------------------------------------------------------------------------------

@stable default_mode = "error" begin
    # Nim-product of the bit patterns a, b < 2^(2^K). Level 0 is GF(2): the product is AND.
    @inline nim_mul_bits(a::T, b::T, ::Val{0}) where {T <: Unsigned} = a & b
    @inline function nim_mul_bits(a::T, b::T, ::Val{K}) where {T <: Unsigned, K}
        h = 1 << (K - 1)                      # half the degree
        lo = (one(T) << h) - one(T)
        a1, a0 = a >> h, a & lo
        b1, b0 = b >> h, b & lo
        c = nim_mul_bits(a1, b1, Val(K - 1))
        d = nim_mul_bits(a0, b0, Val(K - 1))
        e = nim_mul_bits(a0 ⊻ a1, b0 ⊻ b1, Val(K - 1))
        return ((e ⊻ d) << h) | (d ⊻ nim_mul_half(c, Val(K - 1)))
    end

    # c ⊗ 2^(2^k - 1) for c < 2^(2^k): the product with the constant F/2 of the level above,
    # cheaper than a general product. With F′ = 2^(2^(k-1)) one has 2^(2^k - 1) = F′ ⊗ (F′/2), and
    # for c = c1 F′ + c0: c ⊗ F′ = (c1 ⊕ c0) F′ ⊕ c1 ⊗ F′/2, so
    #   H_k(c) = H_{k-1}(c1 ⊕ c0) F′ + H_{k-1}(H_{k-1}(c1)),   H_0(c) = c.
    @inline nim_mul_half(c::T, ::Val{0}) where {T <: Unsigned} = c
    @inline function nim_mul_half(c::T, ::Val{k}) where {T <: Unsigned, k}
        h = 1 << (k - 1)
        lo = (one(T) << h) - one(T)
        c1, c0 = c >> h, c & lo
        return (nim_mul_half(c1 ⊻ c0, Val(k - 1)) << h) |
               nim_mul_half(nim_mul_half(c1, Val(k - 1)), Val(k - 1))
    end

    # Nim-inverse of a ≠ 0. For a = a1 F + a0, the conjugate over the subfield is
    # a1 F + (a0 ⊕ a1) and the norm N = a1² ⊗ F/2 ⊕ a0 ⊗ a1 ⊕ a0² lies in the subfield, so
    # a⁻¹ = (a1 ⊗ N⁻¹) F + ((a0 ⊕ a1) ⊗ N⁻¹).
    @inline nim_inv_bits(a::T, ::Val{0}) where {T <: Unsigned} = a
    @inline function nim_inv_bits(a::T, ::Val{K}) where {T <: Unsigned, K}
        h = 1 << (K - 1)
        lo = (one(T) << h) - one(T)
        a1, a0 = a >> h, a & lo
        sq1 = nim_mul_bits(a1, a1, Val(K - 1))
        N = nim_mul_half(sq1, Val(K - 1)) ⊻ nim_mul_bits(a0, a1, Val(K - 1)) ⊻
            nim_mul_bits(a0, a0, Val(K - 1))
        Ninv = nim_inv_bits(N, Val(K - 1))
        return (nim_mul_bits(a1, Ninv, Val(K - 1)) << h) | nim_mul_bits(a0 ⊻ a1, Ninv, Val(K - 1))
    end
end

# --- Tables (built once, at precompile time) ------------------------------------------------

const NIM256_MUL = [nim_mul_bits(UInt8(a), UInt8(b), Val(3)) for a in 0:255, b in 0:255]
const NIM256_INV = [a == 0 ? 0x00 : nim_inv_bits(UInt8(a), Val(3)) for a in 0:255]

# The smallest generator of the multiplicative group of GF(256), and its log/exp tables.
# EXP has 510 entries so that EXP[log a + log b + 1] needs no reduction modulo 255.
const NIM256_GENERATOR = let
    order(g) = (x = g; k = 1; while x != 0x01; x = NIM256_MUL[x + 1, g + 1]; k += 1; end; k)
    UInt8(findfirst(g -> order(UInt8(g)) == 255, 2:255) + 1)
end
const NIM256_EXP = let e = zeros(UInt8, 510), x = 0x01
    for i in 1:510
        e[i] = x
        x = NIM256_MUL[x + 1, NIM256_GENERATOR + 1]
    end
    e
end
const NIM256_LOG = let l = zeros(UInt8, 256)
    for i in 0:254
        l[NIM256_EXP[i + 1] + 1] = UInt8(i)
    end
    l
end

# --- mul and inverse, per strategy ----------------------------------------------------------

"""
    mul(strategy::MulStrategy, a, b)

The product `a * b` computed with `strategy`. `a * b` uses [`default_strategy`](@ref).
"""
function mul end

"""
    inverse(strategy::MulStrategy, a)

The multiplicative inverse of `a` computed with `strategy`; throws `DivideError` for zero.
`inv(a)` uses [`default_strategy`](@ref).
"""
function inverse end

@stable default_mode = "error" begin
    mul(::Recursive, a::N, b::N) where {K, T, N <: Nimber{K, T}} =
        unsafe_frombits(N, nim_mul_bits(bits(a), bits(b), Val(K)))

    mul(::MulTable, a::N, b::N) where {N <: Nimber} =
        unsafe_frombits(N, @inbounds NIM256_MUL[Int(bits(a)) + 1, Int(bits(b)) + 1])

    function mul(::LogTable, a::N, b::N) where {N <: Nimber}
        (iszero(a) || iszero(b)) && return zero(N)
        @inbounds s = Int(NIM256_LOG[Int(bits(a)) + 1]) + Int(NIM256_LOG[Int(bits(b)) + 1])
        return unsafe_frombits(N, @inbounds NIM256_EXP[s + 1])
    end

    function inverse(::Recursive, a::N) where {K, T, N <: Nimber{K, T}}
        iszero(a) && throw(DivideError())
        return unsafe_frombits(N, nim_inv_bits(bits(a), Val(K)))
    end

    function inverse(::MulTable, a::N) where {N <: Nimber}
        iszero(a) && throw(DivideError())
        return unsafe_frombits(N, @inbounds NIM256_INV[Int(bits(a)) + 1])
    end

    function inverse(::LogTable, a::N) where {N <: Nimber}
        iszero(a) && throw(DivideError())
        return unsafe_frombits(N, @inbounds NIM256_EXP[255 - Int(NIM256_LOG[Int(bits(a)) + 1]) + 1])
    end

    Base.:*(a::N, b::N) where {N <: Nimber} = mul(default_strategy(N), a, b)
    Base.inv(a::Nimber) = inverse(default_strategy(typeof(a)), a)
end

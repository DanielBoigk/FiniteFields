# Faster strategies for the large nimber levels GF(2^16) … GF(2^128) (backlog J2a):
# `LogTable16` stops the shared recursion one level later, at GF(2^16), with log/exp tables;
# `CLMul` multiplies in a polynomial basis GF(2)[x]/(f) with carry-less multiplication, behind a
# change of basis from the nimber basis. Theory and measurements: nimber_fast.md beside this file.

"""
    LogTable16()

The recursion of [`Recursive`](@ref) down to ``\\mathbb{F}_{2^{16}}``, then
``a \\otimes b = \\exp(\\log a + \\log b)`` with logarithm (128 KiB) and exponential (256 KiB)
tables for a fixed generator of ``\\mathbb{F}_{2^{16}}^\\times``. One level fewer than
[`LogTable`](@ref) and [`MulTable`](@ref): GF(2^128) needs 27 base products instead of 81.
Since the subfields are nested, the same tables serve every level up to GF(2^16).
"""
struct LogTable16 <: MulStrategy end

"""
    CLMul()

Multiplication in a polynomial basis: for ``n = 16, 32, 64, 128`` the nimber field
``\\mathbb{F}_{2^n}`` is identified with ``\\mathbb{F}_2[x]/(f)`` for a fixed sparse irreducible
``f`` ([`poly_modulus`](@ref)) by sending ``x`` to a root ``\\alpha`` of ``f`` in the nimber field.
A product converts both factors to the polynomial basis (a bit-matrix product, by byte tables),
multiplies them carry-lessly (x86 `PCLMULQDQ` when the CPU has it, else a software loop), reduces
modulo ``f`` with shifts, and converts back. This is done at ``\\mathbb{F}_{2^{64}}`` and
``\\mathbb{F}_{2^{128}}``; below, where the conversions cost more than they save, `CLMul`
computes as [`MulTable`](@ref), and likewise on CPUs without `PCLMULQDQ` ([`has_pclmul`](@ref)),
where the software carry-less product is slower than the tables. Inversion uses the norm recursion of [`Recursive`](@ref) with
`CLMul` products.
"""
struct CLMul <: MulStrategy end

# --- LogTable16 -----------------------------------------------------------------------------

# The levels the GF(2^16) tables cover.
const Level16 = Union{Val{0}, Val{1}, Val{2}, Val{3}, Val{4}}

# The smallest generator of GF(2^16)^×: order 65535 = 3 · 5 · 17 · 257.
const NIM16_GENERATOR = let
    function powmod(x::UInt16, e::Int)
        r, b = 0x0001, x
        while e > 0
            isodd(e) && (r = nim_mul_bits(MulTable(), r, b, Val(4)))
            b = nim_mul_bits(MulTable(), b, b, Val(4))
            e >>= 1
        end
        return r
    end
    UInt16(findfirst(g -> all(p -> powmod(UInt16(g), 65535 ÷ p) != 0x0001, (3, 5, 17, 257)),
                     2:65535) + 1)
end
# EXP has 2 · 65535 entries so that EXP[log a + log b + 1] needs no reduction modulo 65535.
const NIM16_EXP = let e = zeros(UInt16, 2 * 65535), x = 0x0001
    for i in eachindex(e)
        e[i] = x
        x = nim_mul_bits(MulTable(), x, NIM16_GENERATOR, Val(4))
    end
    e
end
const NIM16_LOG = let l = zeros(UInt16, 65536)
    for i in 0:65534
        l[NIM16_EXP[i + 1] + 1] = UInt16(i)
    end
    l
end
# log 2^(2^k - 1), k = 0..4: the constants of the half products at the levels the tables cover.
const NIM16_LOG_HALF = ntuple(k -> NIM16_LOG[(1 << ((1 << (k - 1)) - 1)) + 1], 5)

@stable default_mode = "error" begin
    @inline function nim_mul_bits(::LogTable16, a::T, b::T, ::Level16) where {T <: Unsigned}
        (iszero(a) || iszero(b)) && return zero(T)
        @inbounds l = Int(NIM16_LOG[Int(a % UInt16) + 1]) + Int(NIM16_LOG[Int(b % UInt16) + 1])
        return (@inbounds NIM16_EXP[l + 1]) % T
    end
    @inline function nim_mul_half(::LogTable16, c::T, v::Level16) where {T <: Unsigned}
        iszero(c) && return zero(T)
        @inbounds l = Int(NIM16_LOG[Int(c % UInt16) + 1]) + Int(NIM16_LOG_HALF[level(v) + 1])
        return (@inbounds NIM16_EXP[l + 1]) % T
    end
    @inline nim_inv_bits(::LogTable16, a::T, ::Level16) where {T <: Unsigned} =
        (@inbounds NIM16_EXP[65535 - Int(NIM16_LOG[Int(a % UInt16) + 1]) + 1]) % T
end

# inv uses LogTable16 where it measured faster than the default of `*` (nimber_fast.md).
default_inv_strategy(::Type{<:Union{GF2_16, GF2_32}}) = LogTable16()

# --- Carry-less multiplication --------------------------------------------------------------

const PCLMUL = Ref(false)

"""
    has_pclmul()

Whether this CPU has the x86 carry-less multiplication instruction `PCLMULQDQ`, which
[`CLMul`](@ref) then uses; otherwise it falls back to [`clmul64_sw`](@ref). Detected at load time.
"""
has_pclmul() = PCLMUL[]

function detect_pclmul()
    @static if Sys.ARCH === :x86_64
        CPUID = Base.BinaryPlatforms.CPUID
        return CPUID.test_cpu_feature(CPUID.JL_X86_pclmul)
    else
        return false
    end
end

"""
    clmul64_sw(a::UInt64, b::UInt64) -> UInt128

The carry-less product of `a` and `b`, i.e. their product as polynomials over ``\\mathbb{F}_2``,
by a loop over the bits of `b` (the portable fallback of [`clmul64`](@ref)).
"""
@inline function clmul64_sw(a::UInt64, b::UInt64)
    A = UInt128(a)
    r = zero(UInt128)
    for i in 0:63
        r ⊻= ifelse(isodd(b >> i), A << i, zero(UInt128))
    end
    return r
end

@static if Sys.ARCH === :x86_64
    const V2 = NTuple{2, VecElement{UInt64}}
    # The function carries "+pclmul" itself, so it compiles whatever the host target says; it is
    # only executed when `has_pclmul()` (otherwise it would be an illegal instruction).
    @inline function clmul64_hw(a::UInt64, b::UInt64)
        r = Base.llvmcall(("""
            declare <2 x i64> @llvm.x86.pclmulqdq(<2 x i64>, <2 x i64>, i8)
            define <2 x i64> @entry(<2 x i64> %a, <2 x i64> %b) #0 {
                %r = call <2 x i64> @llvm.x86.pclmulqdq(<2 x i64> %a, <2 x i64> %b, i8 0)
                ret <2 x i64> %r
            }
            attributes #0 = { alwaysinline "target-features"="+pclmul" }
            """, "entry"), V2, Tuple{V2, V2},
            (VecElement(a), VecElement(zero(UInt64))), (VecElement(b), VecElement(zero(UInt64))))
        return UInt128(r[1].value) | (UInt128(r[2].value) << 64)
    end
    @inline clmul64(a::UInt64, b::UInt64) = PCLMUL[] ? clmul64_hw(a, b) : clmul64_sw(a, b)
else
    @inline clmul64(a::UInt64, b::UInt64) = clmul64_sw(a, b)
end

@doc """
    clmul64(a::UInt64, b::UInt64) -> UInt128

The carry-less product of `a` and `b`: `PCLMULQDQ` if [`has_pclmul`](@ref), else
[`clmul64_sw`](@ref).
""" clmul64

# --- Polynomial basis -----------------------------------------------------------------------

"""
    poly_modulus(T) -> Tuple

The exponents below ``n = 8\\,\\mathrm{sizeof}(T)`` of the irreducible polynomial
``f = x^n + \\sum_e x^e`` over ``\\mathbb{F}_2`` behind [`CLMul`](@ref) for the nimber field
stored in `T`: ``x^{16} + x^5 + x^3 + x + 1``, ``x^{32} + x^7 + x^3 + x^2 + 1``,
``x^{64} + x^4 + x^3 + x + 1``, ``x^{128} + x^7 + x^2 + x + 1`` (the GHASH polynomial).
Irreducibility is checked by the tests (Rabin's test) and by the construction (the powers of the
root are linearly independent).
"""
poly_modulus(::Type{UInt16}) = (5, 3, 1, 0)
poly_modulus(::Type{UInt32}) = (7, 3, 2, 0)
poly_modulus(::Type{UInt64}) = (4, 3, 1, 0)
poly_modulus(::Type{UInt128}) = (7, 2, 1, 0)

# The change of basis for one level: a root α of f in the nimber field, and byte tables for the
# two linear maps (`to`: nimber bits ↦ coordinates in 1, α, …, α^(n-1); `from`: back).
# tab[b + 1, j + 1] is the image of the byte b in position j.
struct PolyBasis{T <: Unsigned}
    root::T
    to::Matrix{T}
    from::Matrix{T}
end

# Root finding at build time: f splits into linear factors over GF(2^n) (it is irreducible of
# degree n), and Cantor–Zassenhaus with the trace splits it: with Tr(Y) = Σ_{i<n} Y^(2^i), the
# factor gcd(g, Tr(βX) mod g) is the product of X - r over the roots r with Tr(βr) = 0, so it is
# proper unless Tr(βr) takes the same value (0 or 1) at every root r; try random β. Follow the
# smaller factor until it is linear. Polynomials are coefficient vectors (lowest first) over F.
function find_root(::Type{F}, exps) where {F <: Nimber}
    n = degree(F)
    trim!(p) = (while !isempty(p) && iszero(p[end]); pop!(p); end; p)
    mt(x, y) = mul(MulTable(), x, y)          # tables only: CLMul needs what is built here
    monic(p) = (c = inverse(MulTable(), p[end]); [mt(x, c) for x in p])
    function pmod(a, g)                          # g monic
        a = trim!(copy(a))
        dg = length(g) - 1
        for i in length(a):-1:(dg + 1)
            c = a[i]
            iszero(c) && continue
            for j in 0:dg
                a[i - dg + j] += mt(c, g[j + 1])
            end
        end
        return trim!(a[1:min(end, dg)])
    end
    function pdiv(a, g)                          # exact quotient, g monic
        a = copy(a)
        dg = length(g) - 1
        q = zeros(F, length(a) - dg)
        for i in length(a):-1:(dg + 1)
            c = a[i]
            q[i - dg] = c
            for j in 0:dg
                a[i - dg + j] += mt(c, g[j + 1])
            end
        end
        return trim!(q)
    end
    function pgcd(a, b)
        a, b = trim!(copy(a)), trim!(copy(b))
        while !isempty(b)
            b = monic(b)
            a, b = b, pmod(a, b)
        end
        return monic(a)
    end
    padd(a, b) = (length(a) < length(b) && ((a, b) = (b, a)); c = copy(a);
                  for i in eachindex(b); c[i] += b[i]; end; trim!(c))
    function psqmod(a, g)
        s = zeros(F, max(2 * length(a) - 1, 0))
        for i in eachindex(a)
            s[2i - 1] = mt(a[i], a[i])
        end
        return pmod(s, g)
    end
    g = zeros(F, n + 1)
    g[n + 1] = one(F)
    for e in exps
        g[e + 1] = one(F)
    end
    rng = Random.Xoshiro(0x4a32)
    while length(g) > 2
        β = rand(rng, F)
        iszero(β) && continue
        t = pmod([zero(F), β], g)
        s = t
        for _ in 1:(n - 1)
            t = psqmod(t, g)
            s = padd(s, t)
        end
        isempty(s) && continue
        h = pgcd(g, s)
        d = length(h) - 1
        (d == 0 || d == length(g) - 1) && continue
        g = 2d <= length(g) - 1 ? h : monic(pdiv(g, h))
    end
    return bits(g[1])                            # g = X + g0, root g0 (characteristic 2)
end

# Byte tables of the linear map with images `img[k + 1]` of the unit vectors e_k:
# tab[b + 1, j + 1] is the image of the byte b in position j. Each entry is one XOR away from
# the entry of b with its lowest set bit cleared, so a table costs 255 XORs per byte position.
byte_tables(img::AbstractVector{T}) where {T <: Unsigned} =
    byte_tables!(Matrix{T}(undef, 256, sizeof(T)), img)
function byte_tables!(tab::Matrix{T}, img::AbstractVector{T}) where {T <: Unsigned}
    @inbounds for j in 0:(sizeof(T) - 1)
        tab[1, j + 1] = zero(T)
        for b in 1:255
            tab[b + 1, j + 1] = tab[(b & (b - 1)) + 1, j + 1] ⊻ img[8j + trailing_zeros(b) + 1]
        end
    end
    return tab
end

function PolyBasis(::Type{F}) where {F <: Nimber}
    T = storagetype(F)
    n = degree(F)
    α = F(find_root(F, poly_modulus(T)))
    cols = Vector{T}(undef, n)                   # α^i in nimber bits
    p = one(F)
    for i in 1:n
        cols[i] = bits(p)
        p = mul(MulTable(), p, α)
    end
    # Invert the matrix with columns α^i: eliminate to pivots (distinct top bits), tracking the
    # combination of columns; then reduce each unit vector through the pivots, highest first.
    pivots = Tuple{T, T}[]                        # (vector, combination of the α^i)
    for i in 1:n
        v, c = cols[i], one(T) << (i - 1)
        for (pv, pc) in pivots
            if !iszero(v & (one(T) << (8 * sizeof(T) - 1 - leading_zeros(pv))))
                v, c = v ⊻ pv, c ⊻ pc
            end
        end
        iszero(v) && error("powers of the root are dependent: f is not irreducible ($T)")
        push!(pivots, (v, c))
        sort!(pivots; by = x -> leading_zeros(x[1]))  # highest top bit first
    end
    inv_img = map(0:(n - 1)) do k
        v, c = one(T) << k, zero(T)
        for (pv, pc) in pivots
            if !iszero(v & (one(T) << (8 * sizeof(T) - 1 - leading_zeros(pv))))
                v, c = v ⊻ pv, c ⊻ pc
            end
        end
        iszero(v) || error("change of basis: elimination failed ($T)")
        c
    end
    return PolyBasis{T}(bits(α), byte_tables(inv_img), byte_tables(cols))
end

const POLY16 = PolyBasis(GF2_16)
const POLY32 = PolyBasis(GF2_32)
const POLY64 = PolyBasis(GF2_64)
const POLY128 = PolyBasis(GF2_128)
poly_basis(::Type{UInt16}) = POLY16
poly_basis(::Type{UInt32}) = POLY32
poly_basis(::Type{UInt64}) = POLY64
poly_basis(::Type{UInt128}) = POLY128
# The tables as constants of their own: one load fewer per lookup than a field of the struct.
for (T, B) in ((UInt16, :POLY16), (UInt32, :POLY32), (UInt64, :POLY64), (UInt128, :POLY128))
    to, from = Symbol(B, :_TO), Symbol(B, :_FROM)
    @eval const $to = $B.to
    @eval const $from = $B.from
    @eval to_table(::Type{$T}) = $to
    @eval from_table(::Type{$T}) = $from
end

"""
    poly_root(T)

The bit pattern of the root ``\\alpha`` of [`poly_modulus`](@ref)`(T)` in the nimber field
stored in `T`, which identifies that field with ``\\mathbb{F}_2[x]/(f)``.
"""
poly_root(::Type{T}) where {T <: Unsigned} = poly_basis(T).root

# The levels where CLMul works in the polynomial basis. Below GF(2^64) the conversions cost more
# than they save (nimber_fast.md, "Measurements"), so CLMul uses MulTable there; the bases of
# GF(2^16) and GF(2^32) are built and tested all the same (for bulk work in the polynomial basis).
const PolyLevel = Union{Val{6}, Val{7}}
const CLMulTableLevel = Union{TableLevel, Val{4}, Val{5}}

@stable default_mode = "error" begin
    @inline function apply_bytes(tab::Matrix{T}, x::T) where {T <: Unsigned}
        r = zero(T)
        for j in 1:sizeof(T)
            r ⊻= @inbounds tab[Int((x >> (8 * (j - 1))) % UInt8) + 1, j]
        end
        return r
    end

    # Nimber bits ↦ coordinates in the basis 1, α, …, α^(n-1), and back.
    @inline to_poly(x::T) where {T <: Unsigned} = apply_bytes(to_table(T), x)
    @inline from_poly(p::T) where {T <: Unsigned} = apply_bytes(from_table(T), p)

    # hi · x^n + lo reduced modulo f = x^n + r: hi · r folds into lo, and the part of hi · r at
    # degree ≥ n (below x^(n + max e)) folds once more; that suffices since 2 max e < n.
    @inline function poly_reduce(hi::T, lo::T) where {T <: Unsigned}
        n = 8 * sizeof(T)
        acc, ov = lo, zero(T)
        for e in poly_modulus(T)
            acc ⊻= hi << e
            e > 0 && (ov ⊻= hi >> (n - e))
        end
        for e in poly_modulus(T)
            acc ⊻= ov << e
        end
        return acc
    end

    # Product in F_2[x]/(f), operands and result as coordinate bit patterns.
    @inline function poly_mulmod(a::T, b::T) where {T <: Union{UInt16, UInt32, UInt64}}
        p = clmul64(a % UInt64, b % UInt64)
        n = 8 * sizeof(T)
        return poly_reduce((p >> n) % T, p % T)
    end
    @inline function poly_mulmod(a::UInt128, b::UInt128)
        a1, a0 = (a >> 64) % UInt64, a % UInt64
        b1, b0 = (b >> 64) % UInt64, b % UInt64
        l = clmul64(a0, b0)
        h = clmul64(a1, b1)
        m = clmul64(a0 ⊻ a1, b0 ⊻ b1) ⊻ l ⊻ h       # Karatsuba: three 64-bit products
        return poly_reduce(h ⊻ (m >> 64), l ⊻ (m << 64))
    end

    # Without PCLMULQDQ the software carry-less product makes this slower than MulTable
    # (nimber_fast.md), so CLMul then computes as MulTable; `poly_mulmod` itself still works.
    @inline function nim_mul_bits(::CLMul, a::T, b::T, v::PolyLevel) where {T <: Unsigned}
        has_pclmul() || return mul_without_pclmul(a, b, v)
        return from_poly(poly_mulmod(to_poly(a), to_poly(b)))
    end
    # Kept out of line, so that the fast path stays small.
    @noinline mul_without_pclmul(a::T, b::T, v::PolyLevel) where {T <: Unsigned} =
        nim_mul_bits(MulTable(), a, b, v)
    @inline nim_mul_bits(::CLMul, a::T, b::T, v::CLMulTableLevel) where {T <: Unsigned} =
        nim_mul_bits(MulTable(), a, b, v)
    @inline nim_mul_half(::CLMul, c::T, v::Val) where {T <: Unsigned} =
        nim_mul_half(MulTable(), c, v)
    @inline nim_inv_bits(::CLMul, a::T, v::CLMulTableLevel) where {T <: Unsigned} =
        nim_inv_bits(MulTable(), a, v)
end

function __init__()
    PCLMUL[] = detect_pclmul()
    return nothing
end

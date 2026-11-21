# Research implementation of Conway's On_2 below ω^(ω^ω) (backlog R1): finite subfields as towers
# of radical extensions over a finite nimber field, following DiMuro, "On On_p" (arXiv:1108.0962v3),
# Theorem 3.1 and Lemma 1.1.
#
# A tower has a base nimber field B = GF(2^(2^a)) and generators χ_r for odd prime powers r = u^k,
# sorted in ordinal order (by u, then k: χ_3 < χ_9 < … < χ_5 < χ_25 < … < χ_7 < …). Generator j
# satisfies χ^u = a_j with a_j in the tower below it: a = χ_{u^(k-1)} for k ≥ 2, a = α_u for k = 1.
# An element is a flat vector of base scalars, index Σ e_j · stride_j (stride_1 = 1), e_j < u_j.
# By Lemma 1.1 its ordinal is Σ_monomials (monomial)·scalar read with the top generator most
# significant: the ordinal order is lexicographic from the highest flat index down, and inside a
# scalar by its bits (finite nimbers are integers). Nim-addition is XOR of everything.
#
# Frobenius powers x ↦ x^(2^s) are kept as maps: on Σ X_i φ^i (φ the top generator) they are
# Σ F_s(X_i) · γ_i · φ^π(i) with φ^(i·2^s) = γ_i φ^π(i); maps compose, so a large s costs no more
# than a small one (no repeated squaring).

module On2

using FiniteFieldArithmetic
using FiniteFieldArithmetic: bits

export Gen, Tower, element, gen_element, scalar_element, mul, sq, frob1, compose, frobpow,
       apply, norm_to, trace_to, field_degree, is_upower, ordinal_string, tower_degree

# --- Towers ----------------------------------------------------------------------------------

struct Gen{B}
    r::Int                    # the prime power, χ_r
    u::Int                    # the prime: χ_r generates a degree-u step
    rel::Vector{B}            # χ_r^u, as an element of the tower below (dense)
end

struct Tower{B}
    gens::Vector{Gen{B}}
    u::Vector{Int}            # radix of each level
    s::Vector{Int}            # size of an element of the level below (s[1] = 1)
    S::Int                    # size of an element
    # workspace per level: unreduced product (2u - 1)·s, and a temporary of size s
    acc::Vector{Vector{B}}
    tmp::Vector{Vector{B}}
end

function Tower(gens::Vector{Gen{B}}) where {B}
    u = [g.u for g in gens]
    s = ones(Int, length(gens))
    for j in 2:length(gens)
        s[j] = s[j - 1] * u[j - 1]
    end
    S = isempty(gens) ? 1 : s[end] * u[end]
    for (j, g) in enumerate(gens)
        length(g.rel) == s[j] || error("relation of χ_$(g.r) has the wrong size")
    end
    Tower{B}(gens, u, s, S, [zeros(B, (2u[j] - 1) * s[j]) for j in eachindex(gens)],
             [zeros(B, s[j]) for j in eachindex(gens)])
end

levels(T::Tower) = length(T.gens)
basetype(::Tower{B}) where {B} = B
"""Degree of the tower over GF(2)."""
tower_degree(T::Tower) = degree(basetype(T)) * T.S

element(T::Tower{B}) where {B} = zeros(B, T.S)
scalar_element(T::Tower{B}, c) where {B} = (x = element(T); x[1] = B(c); x)
"""The generator χ_r of the tower as an element (monomial with exponent 1)."""
function gen_element(T::Tower{B}, r::Int) where {B}
    j = findfirst(g -> g.r == r, T.gens)
    j === nothing && error("χ_$r is not in the tower")
    x = element(T)
    x[T.s[j] + 1] = one(B)
    return x
end

# --- Multiplication ----------------------------------------------------------------------------

iszeroblock(x) = all(iszero, x)

# z .= x * y at level lvl (z, x, y of length s[lvl] * u[lvl], or 1 at level 0)
function mul!(z::AbstractVector{B}, x::AbstractVector{B}, y::AbstractVector{B}, T::Tower{B},
              lvl::Int) where {B}
    if lvl == 0
        @inbounds z[1] = x[1] * y[1]
        return z
    end
    u, s = T.u[lvl], T.s[lvl]
    acc, tmp = T.acc[lvl], T.tmp[lvl]
    fill!(acc, zero(B))
    if lvl == 1                       # scalars: no recursion
        @inbounds for i in 1:u
            xi = x[i]
            iszero(xi) && continue
            for j in 1:u
                acc[i + j - 1] += xi * y[j]
            end
        end
        a = T.gens[1].rel[1]
        @inbounds for t in (2u - 1):-1:(u + 1)
            acc[t - u] += a * acc[t]
        end
        copyto!(z, 1, acc, 1, u)
        return z
    end
    if u >= KARATSUBA_MIN
        copyto!(acc, kara(x, y, u, s, T, lvl))
    else
        schoolbook!(acc, x, y, u, s, T, lvl)
    end
    a = T.gens[lvl].rel
    @inbounds for t in (2u - 2):-1:u
        blk = view(acc, (t * s + 1):((t + 1) * s))
        iszeroblock(blk) && continue
        mul!(tmp, blk, a, T, lvl - 1)
        off = (t - u) * s
        for q in 1:s
            acc[off + q] += tmp[q]
        end
    end
    copyto!(z, 1, acc, 1, u * s)
    return z
end

# Unreduced product of two polynomials of n blocks (level lvl - 1 elements of size s) into acc.
function schoolbook!(acc, x, y, n, s, T::Tower{B}, lvl) where {B}
    tmp = T.tmp[lvl]
    @inbounds for i in 0:(n - 1)
        xi = view(x, (i * s + 1):((i + 1) * s))
        iszeroblock(xi) && continue
        for j in 0:(n - 1)
            yj = view(y, (j * s + 1):((j + 1) * s))
            iszeroblock(yj) && continue
            mul!(tmp, xi, yj, T, lvl - 1)
            off = (i + j) * s
            for q in 1:s
                acc[off + q] += tmp[q]
            end
        end
    end
    return acc
end

const KARATSUBA_MIN = 24

# Karatsuba on block polynomials (characteristic 2: the middle term is P1 + P0 + P2).
function kara(x, y, n, s, T::Tower{B}, lvl) where {B}
    if n < KARATSUBA_MIN ÷ 2
        return schoolbook!(zeros(B, (2n - 1) * s), x, y, n, s, T, lvl)
    end
    h = n ÷ 2
    m = n - h
    x0, x1 = view(x, 1:(h * s)), view(x, (h * s + 1):(n * s))
    y0, y1 = view(y, 1:(h * s)), view(y, (h * s + 1):(n * s))
    P0 = kara(x0, y0, h, s, T, lvl)
    P2 = kara(x1, y1, m, s, T, lvl)
    xs = copy(x1); @views xs[1:(h * s)] .+= x0
    ys = copy(y1); @views ys[1:(h * s)] .+= y0
    P1 = kara(xs, ys, m, s, T, lvl)
    @views P1[1:length(P0)] .+= P0
    P1 .+= P2
    acc = zeros(B, (2n - 1) * s)
    @views acc[1:length(P0)] .+= P0
    @views acc[(h * s + 1):(h * s + length(P1))] .+= P1
    @views acc[(2h * s + 1):(2h * s + length(P2))] .+= P2
    return acc
end

"""x * y in the tower."""
mul(T::Tower, x, y) = mul!(element(T), x, y, T, levels(T))
mul_lvl(T::Tower{B}, x, y, lvl) where {B} = mul!(zeros(B, lvl == 0 ? 1 : T.s[lvl] * T.u[lvl]), x, y, T, lvl)
sq(T::Tower, x) = mul(T, x, x)

add(x, y) = x .+ y

# --- Frobenius maps --------------------------------------------------------------------------

"""x ↦ x^(2^s) on level lvl of a tower: π, γ for the top generator, and the map one level down."""
struct Frob{B}
    s::Int
    lvl::Int
    perm::Vector{Int}                 # π(i), 0-based exponents
    gam::Vector{Vector{B}}            # γ_i in the level below
    sub::Union{Frob{B}, Nothing}
end

lvlsize(T::Tower, lvl) = lvl == 0 ? 1 : T.s[lvl] * T.u[lvl]

scalar_frob(c, s) = (for _ in 1:mod(s, degree(typeof(c))); c = c * c; end; c)

function apply(T::Tower{B}, F::Frob{B}, x::AbstractVector{B}) where {B}
    lvl = F.lvl
    lvl == 0 && return [scalar_frob(x[1], F.s)]
    u, s = T.u[lvl], T.s[lvl]
    z = zeros(B, u * s)
    for i in 0:(u - 1)
        xi = view(x, (i * s + 1):((i + 1) * s))
        iszeroblock(xi) && continue
        c = apply(T, F.sub, xi)
        g = F.gam[i + 1]
        t = isone(g[1]) && all(iszero, view(g, 2:length(g))) ? c : mul_lvl(T, c, g, lvl - 1)
        p = F.perm[i + 1]
        @views z[(p * s + 1):((p + 1) * s)] .+= t
    end
    return z
end
apply(T::Tower, F::Frob, x::Vector) = apply(T, F, view(x, 1:length(x)))

"""Frobenius x ↦ x^2 at level lvl."""
function frob1(T::Tower{B}, lvl::Int = levels(T)) where {B}
    lvl == 0 && return Frob{B}(1, 0, Int[], Vector{B}[], nothing)
    u, s = T.u[lvl], T.s[lvl]
    one_ = zeros(B, s); one_[1] = one(B)
    perm = [(2i) % u for i in 0:(u - 1)]
    gam = [2i < u ? one_ : copy(T.gens[lvl].rel) for i in 0:(u - 1)]
    return Frob{B}(1, lvl, perm, gam, frob1(T, lvl - 1))
end

"""The map Frob^(s1 + s2) = Frob^s2 ∘ Frob^s1."""
function compose(T::Tower{B}, F1::Frob{B}, F2::Frob{B}) where {B}
    lvl = F1.lvl
    lvl == 0 && return Frob{B}(F1.s + F2.s, 0, Int[], Vector{B}[], nothing)
    u = T.u[lvl]
    perm = [F2.perm[F1.perm[i + 1] + 1] for i in 0:(u - 1)]
    gam = [mul_lvl(T, apply(T, F2.sub, F1.gam[i + 1]), F2.gam[F1.perm[i + 1] + 1], lvl - 1)
           for i in 0:(u - 1)]
    return Frob{B}(F1.s + F2.s, lvl, perm, gam, compose(T, F1.sub, F2.sub))
end

"""Frob^s, built from Frob^1 by doubling (the cache keeps the maps already built)."""
function frobpow(T::Tower{B}, s::Integer, cache = Dict{Int, Frob{B}}()) where {B}
    haskey(cache, s) && return cache[s]
    s >= 1 || error("s ≥ 1")
    F = if s == 1
        frob1(T)
    elseif iseven(s)
        h = frobpow(T, s ÷ 2, cache)
        compose(T, h, h)
    else
        compose(T, frobpow(T, s - 1, cache), frobpow(T, 1, cache))
    end
    cache[s] = F
    return F
end

# Σ or Π_{i < k} Frob^(f·i)(x) by doubling; op is add or mul.
function frob_fold(T::Tower, x, f::Int, k::Int, op, cache)
    acc, c = copy(x), 1                       # acc = fold over i < c
    for bit in reverse(digits(k; base = 2))[2:end]
        acc = op(acc, apply(T, frobpow(T, f * c, cache), acc))
        c *= 2
        if bit == 1
            acc = op(x, apply(T, frobpow(T, f, cache), acc))
            c += 1
        end
    end
    @assert c == k
    return acc
end

"""Norm from the tower (degree E) to GF(2^f), f | E."""
norm_to(T::Tower, x, f::Int, cache = Dict()) = frob_fold(T, x, f, tower_degree(T) ÷ f, (a, b) -> mul(T, a, b), cache)
"""Trace from the tower (degree E) to GF(2^f), f | E."""
trace_to(T::Tower, x, f::Int, cache = Dict()) = frob_fold(T, x, f, tower_degree(T) ÷ f, add, cache)

function powe(T::Tower, x, e::Integer)
    r = scalar_element(T, 1)
    b = x
    while e > 0
        isodd(e) && (r = mul(T, r, b))
        e >>= 1
        e > 0 && (b = mul(T, b, b))
    end
    return r
end

primefactors(n) = [p for p in 2:n if n % p == 0 && all(p % q != 0 for q in 2:isqrt(p))]

"""Degree of x over GF(2): the smallest d | E with x^(2^d) = x."""
function field_degree(T::Tower, x, cache = Dict())
    d = tower_degree(T)
    changed = true
    while changed
        changed = false
        for p in primefactors(d)
            if apply(T, frobpow(T, d ÷ p, cache), x) == x
                d ÷= p
                changed = true
                break
            end
        end
    end
    return d
end

"""
Whether x ≠ 0 is a u-th power in the field χ_u, for a prime u not dividing E (the tower degree;
its prime factors are below u): x^((2^E - 1)/u) = 1 (DiMuro Lemma 3.6), computed as
N(x)^((2^f - 1)/u) with N the norm to GF(2^f), f = ord_u(2) | E.
"""
function is_upower(T::Tower, x, u::Int, f::Int, cache = Dict())
    E = tower_degree(T)
    E % u == 0 && error("u divides the tower degree")
    E % f == 0 || error("f ∤ E")
    N = norm_to(T, x, f, cache)
    return powe(T, N, (big(2)^f - 1) ÷ u) == scalar_element(T, 1)
end

# --- Ordinals --------------------------------------------------------------------------------

function ordinal_string(T::Tower, x)
    terms = String[]
    for idx in length(x):-1:1
        c = x[idx]
        iszero(c) && continue
        e = idx - 1
        mono = String[]
        for j in levels(T):-1:1
            ej = (e ÷ T.s[j]) % T.u[j]
            ej == 0 && continue
            push!(mono, ej == 1 ? "χ$(T.gens[j].r)" : "χ$(T.gens[j].r)^$ej")
        end
        cs = string(Int(bits(c)))
        push!(terms, isempty(mono) ? cs : (isone(c) ? join(mono, "·") : join(mono, "·") * "·" * cs))
    end
    return isempty(terms) ? "0" : join(terms, " + ")
end

end

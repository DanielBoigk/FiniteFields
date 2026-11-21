# R1: checks of the rank-encoding data (data/rank/n_*.txt) produced by rank_encoding.jl.
#  1. the minimal polynomial of each internal generator g is irreducible of degree n;
#  2. for n = 2^k the rank encoding is the nimber encoding: products agree with `Nimber`;
#  3. each inclusion GF(2^m) ⊂ GF(2^n) is an injective ring homomorphism, and inclusions compose
#     (GF(2^l) ⊂ GF(2^m) ⊂ GF(2^n) equals GF(2^l) ⊂ GF(2^n));
#  4. GF(8) ⊂ GF(64) agrees with gf8_as_ordinals.jl (ranks 2, 4 ↦ 28, 40);
#  5. order: for small towers, all elements of GF(2^n) are enumerated inside the tower and sorted
#     by ordinal; their positions must be their ranks (the encoding is the ordinal rank).
# Run: julia --project=julia/FiniteFieldArithmetic.jl research/nimber_closure/check_rank_data.jl

using FiniteFieldArithmetic, Random, Test, Serialization
using FiniteFieldArithmetic: bits, field_modulus

const DIR = joinpath(@__DIR__, "data", "rank")

function load(n)
    d = Dict{String, String}()
    for line in eachline(joinpath(DIR, "n_$n.txt"))
        k, v = split(line, '\t')
        d[k] = v
    end
    C = parse.(UInt128, split(d["C"], ','))
    emb = Dict(parse(Int, split(k)[2]) => parse.(UInt128, split(v, ',')) for (k, v) in d if startswith(k, "embed"))
    (; n, g = parse(UInt128, d["g"]), low = parse(UInt128, d["low"]), C, emb)
end

# GF(2) linear algebra on columns (UInt128): apply M to v, invert M
apply(M, v) = (r = zero(UInt128); for i in eachindex(M); isodd(v >> (i - 1)) && (r ⊻= M[i]); end; r)
function inverse(M, n)
    # columns M[i] = image of e_i; solve for images of e_j under M^{-1}
    A = [(M[i], UInt128(1) << (i - 1)) for i in 1:n]
    piv = Tuple{UInt128, UInt128}[]
    for (r, c) in A
        for (pr, pc) in piv
            !iszero(r & (UInt128(1) << (127 - leading_zeros(pr)))) && (r ⊻= pr; c ⊻= pc)
        end
        iszero(r) && error("singular")
        push!(piv, (r, c))
        sort!(piv; by = x -> leading_zeros(x[1]))
    end
    solve(v) = (c = zero(UInt128); for (pr, pc) in piv
        !iszero(v & (UInt128(1) << (127 - leading_zeros(pr)))) && (v ⊻= pr; c ⊻= pc); end; c)
    [solve(UInt128(1) << (j - 1)) for j in 1:n]
end

# Multiplication in rank coordinates: rank → polynomial basis of g (C^{-1}), multiply in
# GF(2)[x]/(P) (the package's PolyElement), back (C).
struct RankField
    n::Int
    P::Type
    C::Vector{UInt128}
    Cinv::Vector{UInt128}
end
RankField(D) = RankField(D.n, poly_field(D.n, D.low), D.C, inverse(D.C, D.n))
function rmul(F::RankField, x::UInt128, y::UInt128)
    T = FiniteFieldArithmetic.storagetype(F.P)
    a = F.P(apply(F.Cinv, x) % T)
    b = F.P(apply(F.Cinv, y) % T)
    apply(F.C, UInt128(bits(a * b)))
end
mask(n) = n == 128 ? typemax(UInt128) : (UInt128(1) << n) - 1

ns = sort([parse(Int, m.captures[1]) for f in readdir(DIR) for m in (match(r"^n_(\d+)\.txt$", f),) if m !== nothing])
println("data for $(length(ns)) degrees: ", ns == collect(1:128) ? "1 … 128" : string(ns))
D = Dict(n => load(n) for n in ns)
F = Dict{Int, RankField}()

@testset "1. minimal polynomials irreducible" begin
    for n in ns
        @test (F[n] = RankField(D[n])) isa RankField       # poly_field checks irreducibility (Rabin)
        @test D[n].C[1] == 1                                 # g^0 = 1 has rank 1
    end
end

rng = Xoshiro(7)
@testset "2. n = 2^k: rank encoding = nimbers" begin
    for (k, N) in enumerate((GF2, GF4, GF16, GF256, GF2_16, GF2_32, GF2_64, GF2_128))
        n = degree(N)
        n in ns || continue
        for _ in 1:200
            x, y = rand(rng, UInt128) & mask(n), rand(rng, UInt128) & mask(n)
            @test rmul(F[n], x, y) == UInt128(bits(N(x) * N(y)))
        end
    end
end

@testset "3. inclusions are homomorphisms and compose" begin
    for n in ns, (m, E) in D[n].emb
        m in ns || continue
        for _ in 1:30
            x, y = rand(rng, UInt128) & mask(m), rand(rng, UInt128) & mask(m)
            @test apply(E, rmul(F[m], x, y)) == rmul(F[n], apply(E, x), apply(E, y))
            @test apply(E, x ⊻ y) == apply(E, x) ⊻ apply(E, y)
        end
        @test apply(E, UInt128(1)) == 1
        lim = m >= 12 ? 4096 : 1 << m                    # injective on the first ranks
        @test length(Set(apply(E, UInt128(k)) for k in 0:(lim - 1))) == lim
        for (l, El) in D[m].emb                                    # l | m | n
            haskey(D[n].emb, l) || continue
            @test all(apply(E, El[i]) == D[n].emb[l][i] for i in 1:l)
        end
    end
end

@testset "4. GF(8) in GF(64) as in gf8_as_ordinals.jl" begin
    @test D[6].emb[3] == UInt128[1, 28, 40]
end

# 5. Enumerate GF(2^n) inside its tower for small towers (degree ≤ 20) and compare positions in
#    ordinal order with ranks.
include(joinpath(@__DIR__, "On2.jl"))
using .On2
@testset "5. rank = position in ordinal order (small towers)" begin
    for n in ns
        R = deserialize(joinpath(DIR, "basis_$n.jls"))
        T, basis = R.T, R.basis
        E = tower_degree(T)
        E <= 20 || continue
        w = degree(eltype(element(T)))
        key(x) = foldl((acc, i) -> acc | (UInt64(bits(x[i])) << ((i - 1) * w)), 1:length(x); init = UInt64(0))
        Fn = frobpow(T, n)
        fixed = UInt64[]
        for k in UInt64(0):(UInt64(1) << E - 1)
            x = [eltype(element(T))((k >> ((i - 1) * w)) & ((1 << w) - 1)) for i in 1:T.S]
            On2.apply(T, Fn, x) == x && push!(fixed, k)                     # x ∈ GF(2^n)
        end
        @test length(fixed) == 2^n
        sort!(fixed)                                      # ordinal order = integer order of digits
        elem(r) = (x = element(T); for i in 0:(n - 1); isodd(r >> i) && (x = x .+ basis[i + 1][2]); end; x)
        @test all(key(elem(UInt128(r))) == fixed[r + 1] for r in 0:(2^n - 1))
        println("  n = $n: all $(2^n) elements, tower degree $E, ranks = ordinal positions")
    end
end

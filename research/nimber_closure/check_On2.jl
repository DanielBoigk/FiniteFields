# R1: independent checks of the On2.jl machinery (run before trusting compute_excesses.jl).
#  1. tower multiplication: commutative, associative, distributive on random elements; x^(2^E) = x;
#  2. Frobenius maps (compose / frobpow) against repeated squaring;
#  3. the u-th power test via norms against the direct exponent x^((2^E - 1)/u);
#  4. the excesses do not depend on the base field (GF(256) vs GF(2^16) as base).
# Run: julia --project=julia/FiniteFieldArithmetic.jl research/nimber_closure/check_On2.jl

include(joinpath(@__DIR__, "On2.jl"))
using .On2, FiniteFieldArithmetic, Random, Test

function tower(B, spec)            # spec: [(r, u, [deps], scalar)] in ordinal order
    gens = Gen{B}[]
    for (r, u, deps, c) in spec
        T = Tower(copy(gens))
        rel = element(T)
        for d in deps
            rel .+= gen_element(T, d)
        end
        rel[1] += B(c)
        push!(gens, Gen{B}(r, u, rel))
    end
    Tower(gens)
end
# Relations from the computed excesses (compute_excesses.jl) for the generators used here.
spec7(B) = tower(B, [(3, 3, Int[], 2), (7, 7, [3], 1)])                  # ω^3 = 2, χ7^7 = ω + 1
spec23(B) = tower(B, [(5, 5, Int[], 4), (11, 11, [5], 1), (23, 23, [11], 1)])
spec29(B) = tower(B, [(3, 3, Int[], 2), (7, 7, [3], 1), (29, 29, [7], 4)])   # a level with u ≥ 24: Karatsuba
rnd(T, rng) = [rand(rng, eltype(element(T))) for _ in 1:length(element(T))]
function sqpow(T, x, k)
    for _ in 1:k
        x = sq(T, x)
    end
    x
end

rng = Xoshiro(5)
@testset "tower arithmetic" begin
    for T in (spec7(GF256), spec23(GF16), spec29(GF16))
        for _ in 1:5
            x, y, z = rnd(T, rng), rnd(T, rng), rnd(T, rng)
            @test mul(T, x, y) == mul(T, y, x)
            @test mul(T, mul(T, x, y), z) == mul(T, x, mul(T, y, z))
            @test mul(T, x, y .+ z) == mul(T, x, y) .+ mul(T, x, z)
        end
        x = rnd(T, rng)
        E = tower_degree(T)
        E <= 2500 && @test sqpow(T, x, E) == x                     # x^(2^E) = x: a field of degree E
    end
end

@testset "Frobenius maps against squaring" begin
    T = spec7(GF256)
    cache = Dict()
    for s in (1, 2, 3, 5, 8, 13, 50, 167)
        x = rnd(T, rng)
        @test apply(T, frobpow(T, s, cache), x) == sqpow(T, x, s)
    end
    T = spec23(GF16)
    for s in (1, 7, 64, 333)
        x = rnd(T, rng)
        @test apply(T, frobpow(T, s, Dict()), x) == sqpow(T, x, s)
    end
end

@testset "u-th power test: norm against direct exponent" begin
    T = spec7(GF256)                       # degree 168 = 8·3·7
    E = tower_degree(T)
    for (u, f) in ((11, 10), (13, 12), (43, 14), (29, 28))
        E % f == 0 || continue
        for _ in 1:3
            x = rnd(T, rng)
            all(iszero, x) && continue
            direct = On2.powe(T, x, (big(2)^E - 1) ÷ u) == scalar_element(T, 1)
            @test is_upower(T, x, u, f) == direct
        end
    end
end

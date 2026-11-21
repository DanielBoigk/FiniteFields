# J3: inversion strategies Fermat(), ItohTsujii(), Euclid(), for nimbers and polynomial-basis
# fields, and the default of `inv` for polynomial-basis fields.
#
# Oracle: x · y = 1 with the field's own product (tested against independent references in
# nimber*.jl and conway.jl), and agreement with the reference inverse of each representation
# (`Recursive` for nimbers, Fermat's x^(2^n - 2) by square-and-multiply on BigInt polynomials
# for polynomial bases, as in test/conway.jl).

using FiniteFieldArithmetic
using FiniteFieldArithmetic: bits, inverse, Recursive, Fermat, ItohTsujii, Euclid, InvStrategy,
    default_inv_strategy
using Random, Test

const INV_STRATEGIES = (Fermat(), ItohTsujii(), Euclid())
const NIMBER_LEVELS = (GF2, GF4, GF16, GF256, GF2_16, GF2_32, GF2_64, GF2_128)
const POLY_FIELDS = Any[conway_field(n) for n in (1, 2, 3, 4, 7, 8, 13, 16, 31, 32, 33, 63, 64, 65, 92, 127)]
append!(POLY_FIELDS, [GF256_AES, GF2_128_GHASH])

# BigInt reference for polynomial bases: x^(2^n - 2) mod f.
function ref_inv_poly(x::F) where {F}
    n = degree(F)
    f = (big(1) << n) | big(FiniteFieldArithmetic.field_modulus(F))
    deg(a) = iszero(a) ? -1 : ndigits(a; base = 2) - 1
    pmod(c) = (while deg(c) >= n; c = xor(c, f << (deg(c) - n)); end; c)
    function mulmod(u, v)             # (names distinct from the outer locals: a closure
        acc = big(0)                  # assigning `r` would share the outer `r`)
        while !iszero(v)
            isodd(v) && (acc = xor(acc, u))
            u, v = u << 1, v >> 1
        end
        return pmod(acc)
    end
    r, a, e = big(1), big(bits(x)), big(2)^n - 2
    while e > 0
        isodd(e) && (r = mulmod(r, a))
        a = mulmod(a, a)
        e >>= 1
    end
    return r
end

inv_samples(F, rng) = degree(F) <= 10 ? [x for x in elements(F) if !iszero(x)] :
    [x for x in vcat(rand(rng, F, 400), [one(F), F(max_bits_of(F))]) if !iszero(x)]
max_bits_of(F) = (big(1) << degree(F)) - 1

@testset "strategies on nimber fields" begin
    rng = Xoshiro(31)
    for F in NIMBER_LEVELS, s in INV_STRATEGIES
        @testset "$F, $s" begin
            for x in inv_samples(F, rng)
                y = inverse(s, x)
                @test y isa F
                @test y == inverse(Recursive(), x)
            end
            @test_throws DivideError inverse(s, zero(F))
        end
    end
    # all of GF(2^16)
    for s in INV_STRATEGIES
        @test all(inverse(s, x) == inverse(Recursive(), x) for x in elements(GF2_16) if !iszero(x))
    end
end

@testset "strategies on polynomial-basis fields" begin
    rng = Xoshiro(32)
    for F in POLY_FIELDS, s in INV_STRATEGIES
        @testset "$F, $s" begin
            xs = inv_samples(F, rng)
            for x in xs
                y = inverse(s, x)
                @test y isa F
                @test isone(x * y)
            end
            for x in Iterators.take(xs, 60)
                @test big(bits(inverse(s, x))) == ref_inv_poly(x)
            end
            @test_throws DivideError inverse(s, zero(F))
        end
    end
    # all of GF(2^16), Conway
    C = conway_field(16)
    for s in INV_STRATEGIES
        @test all(isone(x * inverse(s, x)) for x in elements(C) if !iszero(x))
    end
end

@testset "inv and / use default_inv_strategy" begin
    rng = Xoshiro(33)
    for F in vcat(collect(NIMBER_LEVELS), POLY_FIELDS)
        s = default_inv_strategy(F)
        @test s isa Union{InvStrategy, FiniteFieldArithmetic.MulStrategy}
        x, y = rand(rng, F), rand(rng, F)
        iszero(x) && (x = one(F))
        @test inv(x) == inverse(s, x)
        @test y / x == y * inverse(s, x)
        @test_throws DivideError inv(zero(F))
        @test (@inferred inv(x)) isa F
        for s2 in INV_STRATEGIES
            @test (@inferred inverse(s2, x)) isa F
        end
    end
end

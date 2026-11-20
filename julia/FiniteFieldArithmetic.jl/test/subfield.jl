# J9: the smallest field of an element (minimal_level, minimal_degree, reduce_field).
#
# Oracle, independent of the implementation: [F_2(x) : F_2] is the smallest d ≥ 1 with
# x^(2^d) = x (the Frobenius fixes exactly F_{2^d}), computed by repeated squaring in the field.
# For nimbers that degree is a power of two, 2^j, and F_2(x) is the nimber level j.

using FiniteFieldArithmetic
using FiniteFieldArithmetic: bits
using Random
using Test

const NIMBER_FIELDS = (GF2, GF4, GF16, GF256, GF2_16, GF2_32, GF2_64, GF2_128)
const LEVEL = Dict(GF2 => 0, GF4 => 1, GF16 => 2, GF256 => 3, GF2_16 => 4, GF2_32 => 5,
                   GF2_64 => 6, GF2_128 => 7)

# Smallest d ≥ 1 with x^(2^d) = x, by squaring.
function frobenius_degree(x)
    y, d = x * x, 1
    while y != x
        y, d = y * y, d + 1
    end
    return d
end

# Elements to check: all of GF(2) … GF(2^16), random ones above, plus the edges of every level
# (2^(2^j) - 1 and 2^(2^j), the largest element of level j and the smallest outside it).
function samples(F)
    degree(F) <= 16 && return collect(elements(F))
    rng = Xoshiro(9)
    xs = [rand(rng, F) for _ in 1:2000]
    # random elements of every subfield, so that small degrees are exercised too
    for G in NIMBER_FIELDS
        degree(G) < degree(F) || continue
        append!(xs, F(rand(rng, G)) for _ in 1:200)
    end
    for j in 0:LEVEL[F]
        top = (big(1) << (1 << j)) - 1
        push!(xs, F(top))
        j < LEVEL[F] && push!(xs, F(top + 1))
    end
    return xs
end

@testset "minimal_degree is the degree of F_2(x) (Frobenius oracle)" begin
    for F in NIMBER_FIELDS
        @testset "$F" begin
            for x in samples(F)
                d = frobenius_degree(x)
                @test minimal_degree(x) == d
                @test minimal_degree(x) == 1 << minimal_level(x)
            end
        end
    end
end

@testset "minimal_level from the value" begin
    @test minimal_level(GF2(0)) == 0
    @test minimal_level(GF2(1)) == 0
    @test minimal_level(GF256(0)) == 0
    @test minimal_level(GF256(1)) == 0
    @test minimal_level(GF256(2)) == 1
    @test minimal_level(GF256(3)) == 1
    @test minimal_level(GF256(4)) == 2
    @test minimal_level(GF256(0x0f)) == 2
    @test minimal_level(GF256(0x10)) == 3
    @test minimal_level(GF2_128(typemax(UInt64))) == 6
    @test minimal_level(GF2_128(UInt128(1) << 64)) == 7
    @test minimal_level(GF2_128(typemax(UInt128))) == 7
    # the level does not depend on the field the value is read in (the levels are nested)
    for F in NIMBER_FIELDS
        @test minimal_level(F(1)) == 0
        @test minimal_degree(F(1)) == 1
        @test minimal_level(F(0)) == 0
    end
    @test minimal_level(GF256(0x10)) isa Int
    @test minimal_degree(GF256(0x10)) isa Int
end

@testset "reduce_field" begin
    @test reduce_field(GF2_128(3)) === GF4(3)
    @test reduce_field(GF2_128(0)) === GF2(0)
    @test reduce_field(GF2_128(1)) === GF2(1)
    @test reduce_field(GF2_64(0xab)) === GF256(0xab)
    @test reduce_field(GF2_64(0x1ab)) === GF2_16(0x1ab)
    @test reduce_field(GF2_128(UInt128(1) << 64)) === GF2_128(UInt128(1) << 64)
    @test reduce_field(GF16(5)) === GF16(5)
    for F in NIMBER_FIELDS, x in Iterators.take(samples(F), 3000)
        r = reduce_field(x)
        @test r == x                                        # the same element, embedded
        @test bits(r) == bits(x)
        @test degree(typeof(r)) == minimal_degree(x)        # in the smallest field
        @test reduce_field(r) === r                         # idempotent
        @test minimal_level(r) == minimal_level(x)
    end
    # arithmetic in the reduced field agrees with arithmetic in the original one
    rng = Xoshiro(10)
    for _ in 1:500
        x, y = GF2_64(rand(rng, GF2_16)), GF2_64(rand(rng, GF256))
        @test reduce_field(x) * reduce_field(y) == x * y
        iszero(y) || @test reduce_field(x) / reduce_field(y) == x / y
    end
end

@testset "type stability where promised" begin
    for F in NIMBER_FIELDS
        x = F(1)
        @test (@inferred minimal_level(x)) == 0
        @test (@inferred minimal_degree(x)) == 1
    end
    # reduce_field's result type depends on the value (maintainer, 2026-11-20): not inferable
    @test !isconcretetype(only(Base.return_types(reduce_field, (GF2_64,))))
end

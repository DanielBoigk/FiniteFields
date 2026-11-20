# The default strategy of `inv` (Status "Next" item 1, 2026-11-20): `inv` has its own default,
# `default_inv_strategy`, separate from `default_strategy` behind `*`. Checked here:
# - the defaults per level (LogTable16 at GF(2^16) and GF(2^32), as `*` elsewhere);
# - `inv` and `/` agree with the reference `Recursive` (exhaustively on GF(2^16), on random
#   elements above), zero throws, and the results are type stable.
using FiniteFieldArithmetic
using FiniteFieldArithmetic: bits, inverse, default_strategy, default_inv_strategy, Recursive,
    MulTable, LogTable16, CLMul
using Random, Test

const INV_LEVELS = (GF2, GF4, GF16, GF256, GF2_16, GF2_32, GF2_64, GF2_128)

@testset "defaults" begin
    @test default_inv_strategy(GF2_16) === LogTable16()
    @test default_inv_strategy(GF2_32) === LogTable16()
    for F in (GF2, GF4, GF16, GF256, GF2_64, GF2_128)
        @test default_inv_strategy(F) === default_strategy(F)
    end
    @test default_inv_strategy(GF2_128) === CLMul()
    @test default_strategy(GF2_16) === MulTable()        # `*` is unchanged
end

@testset "inv agrees with Recursive" begin
    @test all(inv(x) == inverse(Recursive(), x) for x in elements(GF2_16) if !iszero(x))
    rng = Xoshiro(2028)
    for F in INV_LEVELS
        xs = [x for x in rand(rng, F, 2000) if !iszero(x)]
        @test all(inv(x) == inverse(Recursive(), x) for x in xs)
        @test all(isone(x * inv(x)) for x in xs)
        y = rand(rng, F)
        @test all(y / x == y * inverse(Recursive(), x) for x in xs)
        @test inv.(xs) == inverse.(Recursive(), xs)
        @test_throws DivideError inv(zero(F))
        @test (@inferred inv(one(F))) isa F
        @test (@inferred default_inv_strategy(F)) isa FiniteFieldArithmetic.MulStrategy
    end
end

using FiniteFields
using Test

@testset "Prime-field arithmetic" begin
    F = Fp{7}
    @test F(8) == F(1)
    @test F(-1) == F(6)
    @test F(3) + F(5) == F(1)
    @test F(3) - F(5) == F(5)
    @test F(3) * F(5) == F(1)
    @test inv(F(3)) == F(5)
    @test F(3) / F(5) == F(2)
    @test zero(F) + one(F) == one(F)
end

@testset "Field validity" begin
    @test_throws ArgumentError Fp{1}(0)
    @test_throws ArgumentError Fp{4}(0)
    @test_throws DivideError inv(Fp{5}(0))
end

# JET.jl: no possible errors found by abstract interpretation in the package, and no runtime
# dispatch in the arithmetic hot paths. (DispatchDoctor's `@stable` in the source additionally
# throws on any type-unstable call to an arithmetic method; the test suite exercises them all.)
using FiniteFieldArithmetic
using FiniteFieldArithmetic: mul, inverse, Recursive, MulTable, LogTable
using JET, Test

@testset "JET" begin
    JET.test_package(FiniteFieldArithmetic; target_modules = (FiniteFieldArithmetic,))
    a, b = GF256(0x53), GF256(0xca)
    for s in (Recursive(), MulTable(), LogTable())
        @test_opt mul(s, a, b)
        @test_opt inverse(s, a)
    end
    @test_opt a + b
    @test_opt a * b
    @test_opt a / b
    @test_opt a^7
    @test_opt GF16(3) * GF256(5)
end

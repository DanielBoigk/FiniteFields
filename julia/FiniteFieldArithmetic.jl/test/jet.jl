# JET.jl: no possible errors found by abstract interpretation in the package, and no runtime
# dispatch in the arithmetic hot paths. (DispatchDoctor's `@stable` in the source additionally
# throws on any type-unstable call to an arithmetic method; the test suite exercises them all.)
using FiniteFieldArithmetic
using FiniteFieldArithmetic: mul, inverse, Recursive, MulTable, LogTable, LogTable16, CLMul,
    ScalarMatMul, TableMatMul, CLMulMatMul, matmul!, mul_elementwise!, ScalarConstMul,
    LazyMulBy, const_mul_broadcasted
using JET, Test

@testset "JET" begin
    JET.test_package(FiniteFieldArithmetic; target_modules = (FiniteFieldArithmetic,))
    a, b = GF256(0x53), GF256(0xca)
    for s in (Recursive(), MulTable(), LogTable(), LogTable16(), CLMul())
        @test_opt mul(s, a, b)
        @test_opt inverse(s, a)
    end
    @test_opt a + b
    @test_opt a * b
    @test_opt a / b
    @test_opt a^7
    @test_opt GF16(3) * GF256(5)
    for F in (GF2_16, GF2_32, GF2_64, GF2_128)
        x, y = F(0x53), F(0xca)
        for s in (Recursive(), MulTable(), LogTable(), LogTable16(), CLMul())
            @test_opt mul(s, x, y)
            @test_opt inverse(s, x)
        end
        @test_opt x^7
    end
    # bulk operations (J8)
    for F in (GF256, GF2_16, GF2_32, GF2_64, GF2_128)
        x = F(0x53)
        A, B, C = fill(x, 3, 3), fill(x, 3, 3), fill(x, 3, 3)
        m = MulBy(x)
        @test_opt m(x)
        for s in (ScalarMatMul(MulTable()), TableMatMul(), CLMulMatMul())
            @test_opt matmul!(s, C, A, B)
        end
        @test_opt A * B
        # constant times array (J8c)
        a = fill(x, 70)
        for s in (ScalarConstMul(), LazyMulBy(0), LazyMulBy(typemax(Int)))
            @test_opt const_mul_broadcasted(s, x, a)
            @test_opt Base.Broadcast.materialize(const_mul_broadcasted(s, x, a))
        end
        @test_opt m(GF2(1))                                       # fallback for other operands
        @test_opt Base.Broadcast.materialize!(similar(a),
            Base.Broadcast.broadcasted(+, const_mul_broadcasted(LazyMulBy(50), x, a), a))
        @test_opt x * a                                # Base broadcasts: the same hook as x .* a
    end
    for F in (GF256, GF2_16, GF2_32)
        v = BitslicedVector(fill(F(3), 70))
        @test_opt mul_elementwise!(v, v, v)
    end
end

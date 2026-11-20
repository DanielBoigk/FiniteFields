# Bulk and array operations (backlog J8): random sampling of arrays, broadcasting, multiplication
# by a fixed constant (`MulBy`), matrix products over the nimber fields (`matmul!`, `*`, `mul!`),
# and bitsliced vectors (`BitslicedVector`).
#
# Oracles, independent of the bulk code:
# - every product is compared with `mul(Recursive(), …)`, the table-free reference;
# - matrix products with the textbook triple sum of `Recursive` products;
# - bitsliced vectors with the element-wise products of the plain vectors;
# - structural checks: associativity of the matrix product, the identity matrix, linearity of
#   multiplication by a constant.

using FiniteFieldArithmetic
using FiniteFieldArithmetic: bits, mul, Recursive, MulTable, LogTable16, CLMul, has_pclmul,
    MatMulStrategy, ScalarMatMul, TableMatMul, CLMulMatMul, matmul!, default_matmul_strategy,
    mul_elementwise!
using LinearAlgebra, Random, Test

rng = Xoshiro(20261121)
const LEVELS = (GF2, GF4, GF16, GF256, GF2_16, GF2_32, GF2_64, GF2_128)
const MATMUL_STRATEGIES = (ScalarMatMul(Recursive()), ScalarMatMul(MulTable()),
                           ScalarMatMul(LogTable16()), TableMatMul(), CLMulMatMul())

refmul(a, b) = mul(Recursive(), a, b)
function refmatmul(A::AbstractMatrix{F}, B::AbstractVecOrMat{F}) where {F}
    C = zeros(F, size(A, 1), size(B, 2))
    for i in axes(A, 1), j in axes(B, 2), k in axes(A, 2)
        C[i, j] += refmul(A[i, k], B[k, j])
    end
    return B isa AbstractVector ? vec(C) : C
end
# Random elements with the extreme values mixed in.
function samples(F, n)
    x = rand(rng, F, n)
    n >= 3 && (x[1] = zero(F); x[2] = one(F); x[3] = F(FiniteFieldArithmetic.max_bits(F)))
    return x
end

@testset "random arrays, $F" for F in LEVELS
    T = FiniteFieldArithmetic.storagetype(F)
    x = rand(Xoshiro(7), F, 1000)
    @test x isa Vector{F} && length(x) == 1000
    @test x == rand(Xoshiro(7), F, 1000)                     # reproducible
    @test all(v -> bits(v) <= FiniteFieldArithmetic.max_bits(F), x)
    M = rand(rng, F, 3, 4)
    @test M isa Matrix{F} && size(M) == (3, 4)
    v = zeros(F, 20)
    rand!(rng, view(v, 2:2:20))                               # non-contiguous: generic path
    @test all(iszero, v[1:2:20])
    if degree(F) <= 8                                         # every element shows up
        @test length(unique(rand(rng, F, 4000))) == 2^degree(F)
    end
    if degree(F) >= 8                                         # every bit varies
        y = bits.(rand(rng, F, 200))
        @test reduce(|, y) == FiniteFieldArithmetic.max_bits(F)
        @test reduce(&, y) == zero(T)
    end
end

@testset "broadcasting, $F" for F in LEVELS
    a, b = samples(F, 100), samples(F, 100)
    c = rand(rng, F)
    @test a .* b == refmul.(a, b)
    @test c .* a == refmul.(Ref(c), a)
    @test a .+ b == [x + y for (x, y) in zip(a, b)]
    @test a .^ 3 == [refmul(x, refmul(x, x)) for x in a]
    nz = filter(!iszero, a)
    @test all(isone, nz .* inv.(nz))
    @test sum(a) == foldl(+, a; init = zero(F))
    out = similar(a)
    f!(out, a, b, c) = (out .= a .* b .+ c; out)
    f!(out, a, b, c)
    @test (@allocated f!(out, a, b, c)) == 0                  # fused, no temporaries
    @test out == refmul.(a, b) .+ c
    @test isbitstype(F)                                       # usable in GPU arrays (not tested)
end

@testset "multiplication by a constant (MulBy), $F" for F in LEVELS
    for c in (zero(F), one(F), samples(F, 6)...)
        m = MulBy(c)
        x = samples(F, 300)
        @test m.(x) == refmul.(Ref(c), x)
        x1, x2 = rand(rng, F), rand(rng, F)
        @test m(x1 + x2) == m(x1) + m(x2)                     # F_2-linear
    end
end

@testset "matrix products, $F" for F in LEVELS
    @test default_matmul_strategy(F) isa MatMulStrategy
    shapes = ((1, 1, 1), (7, 5, 3), (4, 9, 6), (16, 16, 16), (0, 3, 2), (3, 0, 2), (3, 2, 0))
    for (m, p, q) in shapes
        A, B = samples(F, m * p), samples(F, p * q)
        A, B = reshape(A, m, p), reshape(B, p, q)
        R = refmatmul(A, B)
        @test A * B == R                                       # default strategy via mul!
        for s in MATMUL_STRATEGIES
            C = fill(one(F), m, q)                             # overwritten, not accumulated
            @test matmul!(s, C, A, B) === C
            @test C == R
        end
        x = samples(F, p)
        @test A * x == refmatmul(A, x)
    end
    A, B, C = rand(rng, F, 6, 6), rand(rng, F, 6, 6), rand(rng, F, 6, 6)
    @test (A * B) * C == A * (B * C)
    Id = [i == j ? one(F) : zero(F) for i in 1:6, j in 1:6]
    @test A * Id == A && Id * A == A
    # mul!(C, A, B, α, β) = A B α + C β, with α, β Booleans or field elements.
    C0 = rand(rng, F, 6, 6)
    @test mul!(copy(C0), A, B, true, true) == A * B + C0
    @test mul!(copy(C0), A, B, true, false) == A * B
    @test mul!(copy(C0), A, B, false, true) == C0
    α, β = rand(rng, F), rand(rng, F)
    @test mul!(copy(C0), A, B, α, β) == refmul.(A * B, Ref(α)) + refmul.(C0, Ref(β))
    # views (strided, non-contiguous)
    Av, Bv = view(rand(rng, F, 9, 8), 1:2:9, 2:7), view(rand(rng, F, 6, 12), :, 1:3:12)
    @test Av * Bv == refmatmul(Av, Bv)
    for s in MATMUL_STRATEGIES
        @test matmul!(s, similar(Av, 5, 4), Av, Bv) == refmatmul(Av, Bv)
    end
    @test_throws DimensionMismatch matmul!(TableMatMul(), zeros(F, 2, 2), zeros(F, 2, 3), zeros(F, 2, 2))
end

@testset "CLMulMatMul without PCLMULQDQ" begin
    hw = FiniteFieldArithmetic.PCLMUL[]
    try
        FiniteFieldArithmetic.PCLMUL[] = false
        for F in (GF2_16, GF2_64, GF2_128)
            A, B = rand(rng, F, 5, 7), rand(rng, F, 7, 4)
            @test matmul!(CLMulMatMul(), similar(A, 5, 4), A, B) == refmatmul(A, B)
        end
    finally
        FiniteFieldArithmetic.PCLMUL[] = hw
    end
end

@testset "bitsliced vectors, $F" for F in (GF2, GF4, GF16, GF256, GF2_16, GF2_32)
    for n in (0, 1, 63, 64, 65, 200, 1031)
        a, b = samples(F, n), samples(F, n)
        A, B = BitslicedVector(a), BitslicedVector(b)
        @test A isa AbstractVector{F} && length(A) == n
        @test Vector(A) == a && collect(B) == b               # round trip
        n > 0 && @test A[1] == a[1] && A[end] == a[end]
        P = A .* B
        @test P isa BitslicedVector{F}
        @test Vector(P) == refmul.(a, b)
        @test Vector(A .+ B) == a .+ b
        C = BitslicedVector(zeros(F, n))
        @test mul_elementwise!(C, A, B) === C && Vector(C) == refmul.(a, b)
        if n > 0
            x = rand(rng, F)
            C[n] = x
            @test C[n] == x && Vector(C)[1:(n - 1)] == refmul.(a, b)[1:(n - 1)]
        end
        # the padding planes stay zero, so products of the padding never leak into the elements
        @test Vector(BitslicedVector(a) .* BitslicedVector(b)) == Vector(P)
    end
    @test_throws DimensionMismatch BitslicedVector(zeros(F, 3)) .* BitslicedVector(zeros(F, 4))
    @test mul_elementwise!(similar(samples(F, 5)), fill(one(F), 5), fill(one(F), 5)) == fill(one(F), 5)
end

@testset "BitslicedVector is limited to GF(2^32)" begin
    @test_throws ArgumentError BitslicedVector(zeros(GF2_64, 3))
end

# --- J8b: packing by bit-matrix transposes (src/bulk.md, "Bitsliced vectors") ----------------
#
# Oracle: the transpose read off bit by bit.

using FiniteFieldArithmetic: transpose64!

@testset "64 × 64 bit-matrix transpose" begin
    getbit(x, r, c) = (x[r + 1] >> c) & 1
    for x0 in (rand(rng, UInt64, 64), zeros(UInt64, 64), fill(typemax(UInt64), 64),
               [UInt64(1) << r for r in 0:63], [UInt64(1) << (63 - r) for r in 0:63])
        x = copy(x0)
        @test transpose64!(x) === x
        @test all(getbit(x, r, c) == getbit(x0, c, r) for r in 0:63, c in 0:63)
        @test transpose64!(x) == x0                         # an involution
    end
end

@testset "packing from other vectors, $F" for F in (GF2, GF16, GF256, GF2_16, GF2_32)
    a = samples(F, 300)
    @test Vector(BitslicedVector(view(a, 5:250))) == a[5:250]
    @test Vector(BitslicedVector(view(a, 300:-3:1))) == a[300:-3:1]
    P = BitslicedVector(a)
    @test all(P[i] == a[i] for i in eachindex(a))            # getindex agrees with unpacking
    @test (BitslicedVector(a); @allocated BitslicedVector(a)) <= 2 * sizeof(P.planes) + 1024
    out = similar(a, length(a))
    @test (Vector(P); @allocated Vector(P)) <= 2 * sizeof(out) + 1024
end

# Nimber fields above GF(256): GF(2^16), GF(2^32), GF(2^64), GF(2^128) in UInt16 … UInt128.
#
# Exhaustive checks end at 8 bits, so the oracles here are structural and independent of the
# implementation:
# - subfield nesting: values below 2^(2^j) multiply as in GF(2^(2^j)), down to the exhaustively
#   checked GF(256);
# - the Fermat-power rules at every level: F ⊗ x = F·x for x < F, F ⊗ F = 3F/2;
# - the field axioms on random elements;
# - the norm to the half-size subfield, x ⊗ x^(2^h) (h = half the degree), lies in that subfield;
#   the trace to GF(2), Σ x^(2^i), lies in {0, 1};
# - inverses exhaustively on GF(2^16), on random elements above.

using FiniteFieldArithmetic
using FiniteFieldArithmetic: bits, mul, inverse, Recursive, MulTable, LogTable, LogTable16, CLMul
using Random, Test

const LARGE = ((GF2_16, UInt16, 16), (GF2_32, UInt32, 32), (GF2_64, UInt64, 64),
               (GF2_128, UInt128, 128))
const ALL_STRATEGIES = (Recursive(), MulTable(), LogTable(), LogTable16(), CLMul())
rng = Xoshiro(20261119)

# x^(2^k) by k squarings (independent of `^`).
frob(x, k) = (for _ in 1:k; x = x * x; end; x)

@testset "types" begin
    @test GF2_16 === Nimber{4, UInt16}
    @test GF2_32 === Nimber{5, UInt32}
    @test GF2_64 === Nimber{6, UInt64}
    @test GF2_128 === Nimber{7, UInt128}
    for (F, T, n) in LARGE
        @test isbitstype(F) && sizeof(F) == sizeof(T)      # the smallest type that fits
        @test degree(F) == n && characteristic(F) == 2
        @test bits(rand(rng, F)) isa T
    end
    @test_throws ArgumentError Nimber{8, UInt128}(1)        # GF(2^256): no 256-bit storage yet
    @test_throws ArgumentError Nimber{4, UInt32}(1)         # wrong storage for the level
end

@testset "construction" begin
    @test bits(GF2_16(0xbeef)) === 0xbeef
    @test bits(GF2_128(typemax(UInt128))) === typemax(UInt128)
    @test_throws ArgumentError GF2_16(65536)
    @test_throws ArgumentError GF2_32(-1)
    @test_throws ArgumentError GF2_64(UInt128(1) << 64)
    @test GF2_128(GF256(0x1b)) == GF2_128(0x1b)             # embedding keeps the bits
    @test_throws ArgumentError GF2_16(GF2_32(0x10000))      # restriction checks the range
end

@testset "subfields are nested as integers" begin
    for (F, T, n) in LARGE, _ in 1:2000
        a, b = rand(rng, UInt8), rand(rng, UInt8)
        @test bits(F(a) * F(b)) == bits(GF256(a) * GF256(b))
    end
    for (F, T, n) in LARGE[2:end], _ in 1:2000
        a, b = rand(rng, UInt16), rand(rng, UInt16)
        @test bits(F(a) * F(b)) == bits(GF2_16(a) * GF2_16(b))
    end
    for _ in 1:2000
        a, b = rand(rng, UInt64), rand(rng, UInt64)
        @test bits(GF2_128(a) * GF2_128(b)) == bits(GF2_64(a) * GF2_64(b))
    end
    @test promote_type(GF256, GF2_64) === GF2_64
    @test (@inferred GF16(3) * GF2_128(5)) == GF2_128(3) * GF2_128(5)
    @test GF2_32(7) == GF16(7) && hash(GF2_32(7)) == hash(GF16(7))
end

@testset "Fermat-power rules at every level" begin
    for (F, T, n) in LARGE, j in 0:(trailing_zeros(n) - 1)
        Fj = one(T) << (1 << j)                              # F = 2^(2^j) < 2^n
        for _ in 1:200
            x = rand(rng, T) % Fj
            @test bits(F(Fj) * F(x)) == Fj * x               # ordinary product, no overflow
        end
        @test bits(F(Fj) * F(Fj)) == Fj + Fj >> 1            # 3F/2
    end
end

@testset "strategies agree on random elements" begin
    for (F, T, n) in LARGE, _ in 1:3000
        a, b = rand(rng, F), rand(rng, F)
        r = mul(Recursive(), a, b)
        @test mul(MulTable(), a, b) == r
        @test mul(LogTable(), a, b) == r
        @test mul(LogTable16(), a, b) == r
        @test mul(CLMul(), a, b) == r
        if !iszero(a)
            s = inverse(Recursive(), a)
            @test inverse(MulTable(), a) == s
            @test inverse(LogTable(), a) == s
            @test inverse(LogTable16(), a) == s
            @test inverse(CLMul(), a) == s
        end
    end
end

@testset "field axioms on random elements" begin
    for (F, T, n) in LARGE, _ in 1:5000
        a, b, c = rand(rng, F), rand(rng, F), rand(rng, F)
        @test a * b == b * a
        @test (a * b) * c == a * (b * c)
        @test a * (b + c) == a * b + a * c
        @test a * one(F) == a && iszero(a * zero(F))
        iszero(a) || @test a * inv(a) == one(F)
        iszero(b) || @test a / b * b == a
    end
end

@testset "inverses on all of GF(2^16)" begin
    @test all(x * inv(x) == one(GF2_16) for x in elements(GF2_16) if !iszero(x))
    @test all(inverse(Recursive(), x) == inv(x) for x in elements(GF2_16) if !iszero(x))
end

@testset "Frobenius, norm and trace" begin
    for (F, T, n) in LARGE, _ in 1:300
        x = rand(rng, F)
        @test frob(x, n) == x                                # x^(2^n) = x
        h = n ÷ 2
        N = x * frob(x, h)                                   # norm to GF(2^h)
        @test bits(N) < (one(T) << h)
        tr = x
        y = x
        for _ in 1:(n - 1)
            y = y * y
            tr += y
        end
        @test bits(tr) <= 1                                  # trace to GF(2)
    end
end

@testset "powers" begin
    for (F, T, n) in LARGE, _ in 1:300
        x = rand(rng, F)
        iszero(x) && continue
        @test x^(big(2)^n - 1) == one(F)                     # group order, reduced
        @test x^(big(2)^n - 2) == inv(x)                     # Fermat inverse, via powering
        @test x^-1 == inv(x) && x^-2 == inv(x * x)
        k, m = rand(rng, 0:10^6), rand(rng, 0:10^6)
        @test x^k * x^m == x^(k + m)
        @test x^2 == x * x && x^0 == one(F)
    end
    @test zero(GF2_128)^0 == one(GF2_128)
    @test_throws DivideError zero(GF2_64)^-1
end

@testset "printing" begin
    @test string(GF2_16(0x1b)) == "0x001b"
    @test repr(GF2_32(0xdeadbeef)) == "0xdeadbeef"           # fills its storage: no field label
    @test string(GF2_64(1)) == "0x0000000000000001"
    @test string(GF2_128(typemax(UInt128))) == "0x" * "f"^32
end

@testset "type stability" begin
    for (F, T, n) in LARGE
        a, b = rand(rng, F), rand(rng, F) + one(F)
        for s in ALL_STRATEGIES
            @test (@inferred mul(s, a, b)) isa F
            iszero(b) || @test (@inferred inverse(s, b)) isa F
        end
        @test (@inferred a * b) isa F
        @test (@inferred a^123456789) isa F
    end
end

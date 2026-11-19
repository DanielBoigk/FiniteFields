# The fast strategies for the large nimber levels (backlog J2a): `LogTable16` (GF(2^16) log/exp
# tables as the recursion base) and `CLMul` (carry-less multiplication in a polynomial basis behind
# a change of basis).
#
# Oracles, independent of the implementation:
# - carry-less products against a reference on BigInt bit patterns (schoolbook XOR of shifts),
#   and the hardware instruction against the software fallback;
# - the reduction polynomials f are irreducible over GF(2) (Rabin's test, on BigInt
#   polynomials), and the root α found in the nimber field satisfies f(α) = 0 when evaluated with
#   `Recursive`;
# - the change of basis maps the coordinate vector e_i to α^i (computed with `Recursive`), and
#   the two directions are inverse to each other;
# - products and inverses agree with `Recursive` (exhaustively where feasible).

using FiniteFieldArithmetic
using FiniteFieldArithmetic: bits, mul, inverse, Recursive, MulTable, LogTable, LogTable16, CLMul,
    clmul64_sw, clmul64, has_pclmul, poly_modulus, poly_root, to_poly, from_poly, poly_mulmod,
    NIM16_GENERATOR
using Random, Test

rng = Xoshiro(20261120)
const FAST = (LogTable16(), CLMul())
const POLY_LEVELS = ((GF2_16, UInt16, 16), (GF2_32, UInt32, 32), (GF2_64, UInt64, 64),
                     (GF2_128, UInt128, 128))

# Reference GF(2)[x] arithmetic on BigInt bit patterns.
function ref_clmul(a::Integer, b::Integer)
    a, b, r = big(a), big(b), big(0)
    i = 0
    while b >> i != 0
        isodd(b >> i) && (r ⊻= a << i)
        i += 1
    end
    return r
end
deg(p::BigInt) = ndigits(p; base = 2) - 1
function ref_mod(a::BigInt, f::BigInt)
    d = deg(f)
    while a != 0 && deg(a) >= d
        a ⊻= f << (deg(a) - d)
    end
    return a
end
ref_mulmod(a, b, f) = ref_mod(ref_clmul(a, b), f)
function ref_gcd(a::BigInt, b::BigInt)
    while b != 0
        a, b = b, ref_mod(a, b)
    end
    return a
end
# f as a BigInt: x^n plus the low exponents.
modulus_big(T) = (big(1) << (8 * sizeof(T))) ⊻ reduce(⊻, (big(1) << e for e in poly_modulus(T)))

@testset "carry-less multiplication" begin
    for _ in 1:2000
        a, b = rand(rng, UInt64), rand(rng, UInt64)
        @test clmul64_sw(a, b) == ref_clmul(a, b)
        @test clmul64(a, b) == clmul64_sw(a, b)
    end
    @test clmul64_sw(typemax(UInt64), typemax(UInt64)) == ref_clmul(typemax(UInt64), typemax(UInt64))
    @test clmul64(0x8000000000000000, 0x8000000000000000) == UInt128(1) << 126
    @info "PCLMULQDQ available: $(has_pclmul())"
end

@testset "reduction polynomial and change of basis, GF(2^$n)" for (F, T, n) in POLY_LEVELS
    f = modulus_big(T)
    @test deg(f) == n
    # Rabin: f of degree n = 2^k is irreducible iff x^(2^n) ≡ x mod f and gcd(x^(2^(n/2)) - x, f) = 1.
    sq(p) = ref_mulmod(p, p, f)
    x = big(2)
    y = x
    for _ in 1:(n ÷ 2); y = sq(y); end
    @test ref_gcd(f, y ⊻ x) == 1
    for _ in 1:(n ÷ 2); y = sq(y); end
    @test y == x

    α = F(poly_root(T))
    pow(z, k) = (r = one(F); for _ in 1:k; r = mul(Recursive(), r, z); end; r)
    @test iszero(pow(α, n) + sum(pow(α, e) for e in poly_modulus(T)))       # f(α) = 0
    p = one(F)
    for i in 0:(n - 1)
        @test from_poly(one(T) << i) == bits(p)                             # e_i ↦ α^i
        @test to_poly(bits(p)) == one(T) << i
        p = mul(Recursive(), p, α)
    end
    for _ in 1:1000
        v, w = rand(rng, T), rand(rng, T)
        @test from_poly(to_poly(v)) == v && to_poly(from_poly(v)) == v
        @test to_poly(v ⊻ w) == to_poly(v) ⊻ to_poly(w)
        @test poly_mulmod(v, w) == ref_mulmod(v, w, f)
    end
end

@testset "GF(2^16) generator" begin
    g = GF2_16(NIM16_GENERATOR)
    @test isone(g^65535)
    @test all(!isone(g^(65535 ÷ p)) for p in (3, 5, 17, 257))
end

@testset "agreement with Recursive at the table levels" begin
    for s in FAST, a in elements(GF256), b in elements(GF256)
        @test mul(s, a, b) == mul(MulTable(), a, b)
    end
    for s in FAST, a in elements(GF256)
        iszero(a) || @test inverse(s, a) == inverse(MulTable(), a)
    end
    for s in FAST, F in (GF2, GF4, GF16), a in elements(F), b in elements(F)
        @test mul(s, a, b) == mul(Recursive(), a, b)
    end
end

@testset "GF(2^16): inverses exhaustively, products on a large sample" begin
    for s in FAST
        @test all(mul(s, x, inverse(s, x)) == one(GF2_16) for x in elements(GF2_16) if !iszero(x))
        @test all(inverse(s, x) == inverse(Recursive(), x) for x in elements(GF2_16) if !iszero(x))
        @test all(mul(s, x, y) == mul(Recursive(), x, y)
                  for (x, y) in zip(rand(rng, GF2_16, 200_000), rand(rng, GF2_16, 200_000)))
        @test all(mul(s, x, GF2_16(0x8000)) == mul(Recursive(), x, GF2_16(0x8000))
                  for x in elements(GF2_16))
    end
end

@testset "agreement with Recursive above GF(2^16)" begin
    for (F, T, n) in POLY_LEVELS[2:end], s in FAST, _ in 1:5000
        a, b = rand(rng, F), rand(rng, F)
        @test mul(s, a, b) == mul(Recursive(), a, b)
        iszero(a) || @test inverse(s, a) == inverse(Recursive(), a)
    end
    # Sparse and extreme operands: subfield elements, Fermat powers, all ones.
    for (F, T, n) in POLY_LEVELS, s in FAST
        special = [zero(T), one(T), typemax(T), (one(T) << (n - 1)), (one(T) << (n ÷ 2)),
                   (one(T) << (n ÷ 2)) - one(T), T(0x53), T(0xca)]
        for a in special, b in special
            @test mul(s, F(a), F(b)) == mul(Recursive(), F(a), F(b))
        end
    end
end

@testset "type stability" begin
    for (F, T, n) in POLY_LEVELS, s in FAST
        a, b = rand(rng, F), rand(rng, F) + one(F)
        @test (@inferred mul(s, a, b)) isa F
        iszero(b) || @test (@inferred inverse(s, b)) isa F
    end
end

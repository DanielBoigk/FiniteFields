# Characteristic-2 elements print in hexadecimal (maintainer, 2026-11-19). A field that fills its
# storage shows the bare hex value; one that does not (e.g. GF(16) in a UInt8) also names the
# field, except where the context already states the element type (inside a typed array).

using FiniteFieldArithmetic
using Test

@testset "print is hex with one digit per four bits" begin
    @test string(GF256(0x1b)) == "0x1b"
    @test string(GF256(0x05)) == "0x05"      # padded to the field's width
    @test string(GF16(0xb)) == "0xb"
    @test string(GF4(3)) == "0x3"
    @test string(GF2(1)) == "0x1"
end

@testset "show (REPL)" begin
    @test repr(GF256(0x1b)) == "0x1b"
    @test repr(GF16(0xb)) == "0xb ∈ GF(2^4)"
    @test repr(GF4(2)) == "0x2 ∈ GF(2^2)"
    @test repr(GF2(1)) == "0x1 ∈ GF(2)"
    @test repr("text/plain", GF256(0xff)) == "0xff"
    @test repr("text/plain", GF16(0xb)) == "0xb ∈ GF(2^4)"
end

@testset "arrays do not repeat the field" begin
    @test repr([GF16(1), GF16(0xa)]) == "GF16[0x1, 0xa]"
    @test repr([GF256(1), GF256(0xab)]) == "GF256[0x01, 0xab]"
end

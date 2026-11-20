# Packing and unpacking of BitslicedVectors (J8b, src/bulk.md "Bitsliced vectors"): ns per
# element for 2^14 random elements, three runs per field in one process.
# julia --project=bench bench/pack.jl
using FiniteFieldArithmetic, BenchmarkTools, Random

let rng = Xoshiro(1)
    for F in (GF2, GF16, GF256, GF2_16, GF2_32)
        a = rand(rng, F, 2^14)
        A = BitslicedVector(a)
        for _ in 1:3
            p = @belapsed BitslicedVector($a)
            u = @belapsed Vector($A)
            println(rpad(string(F), 8), " pack ", round(p * 1e9 / 2^14; digits = 3),
                    "  unpack ", round(u * 1e9 / 2^14; digits = 3), " ns per element")
        end
    end
end

# J9: cost of minimal_level / minimal_degree (bit length, type-stable) and of reduce_field
# (type-unstable: the result type depends on the value), on random elements.
# Run: julia --project=julia/FiniteFieldArithmetic.jl/bench julia/FiniteFieldArithmetic.jl/bench/subfield.jl
using FiniteFieldArithmetic, BenchmarkTools, Random

for F in (GF256, GF2_16, GF2_32, GF2_64, GF2_128)
    xs = rand(Xoshiro(1), F, 1024)
    lev = @belapsed sum(minimal_level, $xs)
    red = @belapsed foreach(x -> Base.donotdelete(reduce_field(x)), $xs)
    println(rpad(string(F), 22), "minimal_level ", round(lev / 1024 * 1e9; digits = 2),
        " ns   reduce_field ", round(red / 1024 * 1e9; digits = 1), " ns   (per element)")
end

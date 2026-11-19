# Bulk sampling `rand!(rng, a)` (src/bulk.jl) against one `rand(rng, F)` per element, per element
# over 2^14 elements (src/bulk.md, "Measurements").
#   julia --project=bench bench/rand_compare.jl
using BenchmarkTools, Random, FiniteFieldArithmetic
rng = Xoshiro(1)
for F in (GF16, GF256, GF2_16, GF2_32, GF2_64, GF2_128)
    a = Vector{F}(undef, 2^14)
    t_bulk = @belapsed rand!($rng, $a)
    t_elt = @belapsed (for i in eachindex($a); $a[i] = rand($rng, $F); end)
    println(F, ": rand! ", round(t_bulk / 2^14 * 1e9; digits = 3), " ns, per element ",
            round(t_elt / 2^14 * 1e9; digits = 3), " ns")
end

# Matrix products at growing sizes (J8b, cache blocking; src/bulk.md "Matrix products at larger
# sizes"): ns per multiply-add of an n × n times n × n product with the default strategy,
# minimum of a few runs. About 20 minutes; GF(2^128) at 2048 alone takes about 2.
# julia --project=bench bench/matmul_sizes.jl
using FiniteFieldArithmetic, BenchmarkTools, Random
using FiniteFieldArithmetic: matmul!, default_matmul_strategy

const SIZES = (64, 128, 256, 512, 1024, 2048)
const FIELDS = (GF256, GF2_16, GF2_32, GF2_64, GF2_128)

function measure(F, n)
    rng = Xoshiro(n)
    A, B, C = rand(rng, F, n, n), rand(rng, F, n, n), Matrix{F}(undef, n, n)
    s = default_matmul_strategy(F)
    matmul!(s, C, A, B)                                    # compile, warm caches
    t = @belapsed matmul!($s, $C, $A, $B) samples = (n >= 1024 ? 3 : 20) evals = 1 seconds = 60
    return t * 1e9 / n^3
end

for F in FIELDS, n in SIZES
    println(rpad(string(F), 10), lpad(n, 5), "  ", round(measure(F, n); digits = 3),
            " ns per multiply-add")
    flush(stdout)
end

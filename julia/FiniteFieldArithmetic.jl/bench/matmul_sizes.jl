# Matrix products at growing sizes (J8b, cache blocking; src/bulk.md "Matrix products at larger
# sizes"): ns per multiply-add of an n × n times n × n product with the default strategy,
# minimum of a few runs. About 20 minutes; GF(2^128) at 2048 alone takes about 2.
# julia --project=bench bench/matmul_sizes.jl
# With the argument `compare` (J8b, the bitsliced product; src/bitsliced_matmul.md): every field
# from GF(2) up, each size with `CLMulMatMul` (`TableMatMul` up to GF(256)) and `BitslicedMatMul`,
# the conversions included. About 30 minutes.
using FiniteFieldArithmetic, BenchmarkTools, Random
using FiniteFieldArithmetic: matmul!, default_matmul_strategy, TableMatMul, CLMulMatMul,
    BitslicedMatMul

const SIZES = (64, 128, 256, 512, 1024, 2048)
const FIELDS = (GF256, GF2_16, GF2_32, GF2_64, GF2_128)

function measure(F, n, s = default_matmul_strategy(F))
    rng = Xoshiro(n)
    A, B, C = rand(rng, F, n, n), rand(rng, F, n, n), Matrix{F}(undef, n, n)
    matmul!(s, C, A, B)                                    # compile, warm caches
    t = @belapsed matmul!($s, $C, $A, $B) samples = (n >= 1024 ? 3 : 20) evals = 1 seconds = 60
    return t * 1e9 / n^3
end

if "compare" in ARGS
    for F in (GF2, GF4, GF16, FIELDS...), n in SIZES
        s0 = degree(F) <= 8 ? TableMatMul() : CLMulMatMul()
        t0, t1 = measure(F, n, s0), measure(F, n, BitslicedMatMul())
        println(rpad(string(F), 10), lpad(n, 5), "  ", rpad(string(nameof(typeof(s0))), 12),
                lpad(round(t0; digits = 3), 8), "  BitslicedMatMul ", lpad(round(t1; digits = 3), 8),
                "  ratio ", round(t0 / t1; digits = 2))
        flush(stdout)
    end
else
    for F in FIELDS, n in SIZES
        println(rpad(string(F), 10), lpad(n, 5), "  ", round(measure(F, n); digits = 3),
                " ns per multiply-add")
        flush(stdout)
    end
end

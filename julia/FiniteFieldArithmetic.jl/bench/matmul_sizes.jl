# Matrix products at growing sizes (J8b, cache blocking; src/bulk.md "Matrix products"): ns per
# multiply-add of an n × n times n × n product with the default strategy, minimum of a few runs.
# julia --project=bench bench/matmul_sizes.jl [strategy ...]
#   strategies: default (matmul! with default_matmul_strategy), clmul (CLMulMatMul()),
#   table (TableMatMul()), blocked (CLMulMatMul with blocking, if defined)
using FiniteFieldArithmetic, BenchmarkTools, Random
using FiniteFieldArithmetic: matmul!, default_matmul_strategy

const SIZES = (64, 128, 256, 512, 1024, 2048)
const FIELDS = (GF256, GF2_16, GF2_32, GF2_64, GF2_128)

function measure(F, n, s; maxtime = 60.0)
    rng = Xoshiro(n)
    A, B, C = rand(rng, F, n, n), rand(rng, F, n, n), Matrix{F}(undef, n, n)
    matmul!(s, C, A, B)                                    # compile, warm caches
    # one evaluation per sample; few samples for the large sizes
    t = @belapsed matmul!($s, $C, $A, $B) samples = (n >= 1024 ? 3 : 20) evals = 1 seconds = maxtime
    return t * 1e9 / n^3
end

let strategies = isempty(ARGS) ? ["default"] : ARGS
    for name in strategies, F in FIELDS
        s = name == "default" ? default_matmul_strategy(F) : getfield(FiniteFieldArithmetic, Symbol(name))()
        for n in SIZES
            # skip what takes longer than about 20 s per evaluation
            (F === GF2_128 && n > 1024 && name != "default") && continue
            println(rpad(name, 10), rpad(string(F), 10), lpad(n, 5), "  ",
                    round(measure(F, n, s); digits = 3), " ns per multiply-add")
            flush(stdout)
        end
    end
end

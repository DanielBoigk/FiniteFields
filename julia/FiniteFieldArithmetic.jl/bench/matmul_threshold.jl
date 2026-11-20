# Where BitslicedMatMul starts to win (J8g, src/bitsliced_matmul.md "Default by size"): per field,
# ns per multiply-add of the field's small-size strategy (`TableMatMul` up to GF(256),
# `CLMulMatMul` above) and of `BitslicedMatMul`, conversions included, minimum of several runs.
#
# 1. Square products n × n × n on a finer grid than bench/matmul_sizes.jl; the threshold printed is
#    the smallest n from which BitslicedMatMul was faster at that n and every larger one.
# 2. Rectangular products m × p × q (A is m × p, B is p × q), to check the rule "bitsliced iff
#    min(m, p, q) ≥ threshold": the smallest dimension at the threshold t with the others 4t, and
#    each dimension alone at t/2, 2 and 1 (matrix × vector) with the others at 4t.
#
# julia --project=bench bench/matmul_threshold.jl            (about 30 minutes)
# julia --project=bench bench/matmul_threshold.jl GF2_16     (one field)
using FiniteFieldArithmetic, BenchmarkTools, Random
using FiniteFieldArithmetic: matmul!, TableMatMul, CLMulMatMul, BitslicedMatMul

const GRID = (16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768)
const ALL = (GF2, GF4, GF16, GF256, GF2_16, GF2_32, GF2_64, GF2_128)

small_strategy(F) = degree(F) <= 8 ? TableMatMul() : CLMulMatMul()

function ns_per_muladd(F, s, m, p, q)
    rng = Xoshiro(m + 7p + 13q)
    A, B = rand(rng, F, m, p), rand(rng, F, p, q)
    C = q == 1 ? Vector{F}(undef, m) : Matrix{F}(undef, m, q)
    Bq = q == 1 ? vec(B) : B
    matmul!(s, C, A, Bq)                                   # compile, warm caches
    t = @belapsed matmul!($s, $C, $A, $Bq) samples = 5 evals = 1 seconds = 20
    return t * 1e9 / (m * p * q)
end

function compare(F, m, p, q)
    t0 = ns_per_muladd(F, small_strategy(F), m, p, q)
    t1 = ns_per_muladd(F, BitslicedMatMul(), m, p, q)
    println(rpad(string(F), 9), lpad("$m × $p × $q", 18), "  ", rpad(string(nameof(typeof(small_strategy(F)))), 12),
            lpad(round(t0; sigdigits = 3), 8), "  Bitsliced ", lpad(round(t1; sigdigits = 3), 8),
            "  ratio ", round(t0 / t1; digits = 2))
    flush(stdout)
    return t0, t1
end

function threshold(F)
    wins = [(n, compare(F, n, n, n)...) for n in GRID]
    t = nothing
    for k in length(wins):-1:1
        wins[k][3] < wins[k][2] || break
        t = wins[k][1]
    end
    println(rpad(string(F), 9), "  threshold (square): ", something(t, "none on the grid"))
    return t
end

fields = isempty(ARGS) ? ALL : Tuple(getfield(FiniteFieldArithmetic, Symbol(a)) for a in ARGS)
for F in fields
    t = threshold(F)
    t === nothing && continue
    big = min(4t, 1024)
    for (m, p, q) in ((t, big, big), (big, t, big), (big, big, t),
                      (t ÷ 2, big, big), (big, t ÷ 2, big), (big, big, t ÷ 2),
                      (2, big, big), (big, 2, big), (big, big, 2), (big, big, 1))
        compare(F, m, p, q)
    end
end

# FF-P2: CLMulMatMul for a matrix times a vector or a few columns (src/bulk.md, "Matrix × vector").
# ns per multiply-add, conversions included, minimum of several runs, of the two paths of
# CLMulMatMul:
#
# - `transposed`: A converted to the polynomial basis and transposed once, then one `poly_dot` per
#   entry of the result (the path for many columns);
# - `columns`: column by column, A read in storage order and converted on the fly, unreduced
#   products accumulated per row (the path up to `clmul_column_max(F)` columns);
#
# and, for comparison, `TableMatMul` and `ScalarMatMul(CLMul())`. The last lines per field print
# the largest q of the grid up to which `columns` was faster at every q and size: the measured
# basis of `clmul_column_max`.
#
# julia --project=bench bench/matvec.jl            (about 20 minutes)
# julia --project=bench bench/matvec.jl GF2_64     (one field)
using FiniteFieldArithmetic, BenchmarkTools, Random
using FiniteFieldArithmetic: matmul!, TableMatMul, ScalarMatMul, CLMul, clmul_transposed_add!,
    clmul_columns_add!

const SHAPES = ((1024, 1024), (128, 128), (16, 1024), (1024, 16))
const QS = (1, 2, 3, 4, 6, 8, 12, 16)
const ALL = (GF2_16, GF2_32, GF2_64, GF2_128)

function ns(f!, F, m, p, q)
    rng = Xoshiro(m + 7p + 13q)
    A = rand(rng, F, m, p)
    B = q == 1 ? rand(rng, F, p) : rand(rng, F, p, q)
    C = q == 1 ? Vector{F}(undef, m) : Matrix{F}(undef, m, q)
    f!(C, A, B)                                                  # compile, warm caches
    t = @belapsed $f!(fill!($C, zero($F)), $A, $B) samples = 7 evals = 1 seconds = 3
    return t * 1e9 / (m * p * q)
end

r(x) = lpad(round(x; sigdigits = 3), 7)

function sweep(F)
    best = typemax(Int)
    for (m, p) in SHAPES
        last_win = 0
        for q in QS
            tt = ns(clmul_transposed_add!, F, m, p, q)
            tc = ns(clmul_columns_add!, F, m, p, q)
            extra = ""
            if q == 1 || q == 4
                ttab = ns((C, A, B) -> matmul!(TableMatMul(), C, A, B), F, m, p, q)
                tsc = ns((C, A, B) -> matmul!(ScalarMatMul(CLMul()), C, A, B), F, m, p, q)
                extra = "  Table" * r(ttab) * "  Scalar(CLMul)" * r(tsc)
            end
            println(rpad(string(F), 8), lpad("$m × $p × $q", 18), "  transposed", r(tt),
                    "  columns", r(tc), "  ratio ", lpad(round(tt / tc; digits = 2), 5), extra)
            flush(stdout)
            tc < tt && last_win == (q == 1 ? 0 : QS[findfirst(==(q), QS) - 1]) && (last_win = q)
        end
        println(rpad(string(F), 8), lpad("$m × $p", 18), "  columns faster up to q = ", last_win)
        best = min(best, last_win)
    end
    println(rpad(string(F), 8), "  columns faster at every size up to q = ", best)
    return best
end

fields = isempty(ARGS) ? ALL : Tuple(getfield(FiniteFieldArithmetic, Symbol(a)) for a in ARGS)
foreach(sweep, fields)

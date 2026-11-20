# Multiplication of an array by a constant in broadcasts (backlog J8c, J8d; src/bulk.md): the
# strategies of `c .* a` against each other, per field and array length, alone and fused.
#   julia --project=bench bench/constmul.jl [quick]
#
# Strategies: C = ScalarConstMul(), B = LazyMulBy(0) (tables always), B0 = LazyMulBy(typemax(Int))
# (never tables: the cost of B below its threshold, i.e. of its Union result). The eager
# strategy A of J8c was removed in J8d; its figures are in src/bulk.md.
# Expressions, each written as Julia lowers the dot syntax, with the strategy's hook in place of
# `broadcasted(*, c, a)`:
#   new:   c .* a              materialize(hook(c, a))
#   into:  out .= c .* a       materialize!(out, hook(c, a))
#   fused: out .= c .* a .+ b  materialize!(out, broadcasted(+, hook(c, a), b))
# Output: median ns per element, minimum allocations, and per field the smallest measured length
# from which B beats C in every expression (the break-even that sets `mulby_threshold`).
using BenchmarkTools, Random
using FiniteFieldArithmetic
using FiniteFieldArithmetic: ScalarConstMul, LazyMulBy, const_mul_broadcasted
using Base.Broadcast: broadcasted, materialize, materialize!

const QUICK = "quick" in ARGS
const FIELDS = (GF256, GF2_16, GF2_32, GF2_64, GF2_128)
const LENGTHS = QUICK ? (2^2, 2^6, 2^10, 2^14) :
                (4, 8, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, 2048, 4096,
                 2^14, 2^16)
const STRATS = (("C", ScalarConstMul()), ("B", LazyMulBy(0)), ("B0", LazyMulBy(typemax(Int))))

e_new(s, c, a) = materialize(const_mul_broadcasted(s, c, a))
e_into(s, out, c, a) = materialize!(out, const_mul_broadcasted(s, c, a))
e_fused(s, out, c, a, b) = materialize!(out, broadcasted(+, const_mul_broadcasted(s, c, a), b))

const EXPRS = ("new", "into", "fused")

function measure(F, n, s; seconds = QUICK ? 0.05 : 0.25)
    rng = Xoshiro(n)
    a, b, out, c = rand(rng, F, n), rand(rng, F, n), Vector{F}(undef, n), rand(rng, F)
    r = Dict{String, Tuple{Float64, Int}}()
    for (name, t) in (("new", @benchmarkable e_new($s, $c, $a) seconds = seconds),
                      ("into", @benchmarkable e_into($s, $out, $c, $a) seconds = seconds),
                      ("fused", @benchmarkable e_fused($s, $out, $c, $a, $b) seconds = seconds))
        tune!(t)
        tr = run(t)
        r[name] = (median(tr).time / n, tr.allocs)
    end
    return r
end

println("Julia ", VERSION, ", ", Sys.CPU_NAME, ", ", Sys.CPU_THREADS, " threads")
println("median ns per element (allocations per call)")
for F in FIELDS
    println("\n", F)
    println(rpad("n", 7), join((rpad(sn * " " * e, 19) for (sn, _) in STRATS for e in EXPRS)))
    results = Dict{Tuple{String, Int}, Dict{String, Tuple{Float64, Int}}}()
    for n in LENGTHS
        line = rpad(string(n), 7)
        for (sn, s) in STRATS
            r = results[(sn, n)] = measure(F, n, s)
            for e in EXPRS
                t, al = r[e]
                line *= rpad(string(round(t; digits = 3), " (", al, ")"), 19)
            end
        end
        println(line)
    end
    wins(n) = all(results[("B", n)][e][1] < results[("C", n)][e][1] for e in EXPRS)
    be = findfirst(i -> all(wins, LENGTHS[i:end]), eachindex(LENGTHS))
    println("B beats C in every expression from n = ",
            be === nothing ? "never (in the measured range)" : LENGTHS[be])
end

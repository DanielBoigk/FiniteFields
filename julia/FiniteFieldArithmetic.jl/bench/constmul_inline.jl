# Why the function barrier in LazyMulBy's copyto! (src/bulk.jl, unwrapped_copyto!) pays (backlog
# J8b; src/bulk.md, "Why inlining both branches slowed the loop"). Compares, for the fused
# `out .= c .* a .+ b` below the threshold (Base's `*` branch):
#   C          Base's broadcast (ScalarConstMul)
#   barrier    LazyMulBy as in the package: each branch behind a @noinline function
#   inlined    both branches inlined into one method (J8d's first design)
#   false only only the `*` branch, inlined (no second loop in the method)
# and counts, in the LLVM IR of `inlined` and `false only`, the loop copies and the `select`s left
# in loops by Base's broadcast index (`newindex`: an argument of length 1 is extruded). Run it
# twice, the second time with LLVM's unswitching cost multiplier off:
#   julia --project=bench bench/constmul_inline.jl [field] [n]
#   JULIA_LLVM_ARGS=--enable-unswitch-cost-multiplier=false julia --project=bench bench/constmul_inline.jl
using BenchmarkTools, Random, InteractiveUtils
using FiniteFieldArithmetic
using FiniteFieldArithmetic: ScalarConstMul, LazyMulBy, const_mul_broadcasted, max_threshold,
                             unwrap_constmul
using Base.Broadcast: broadcasted, materialize!, instantiate

function inlined_copyto!(dest, bc)
    return length(dest) >= max_threshold(bc) ?
           copyto!(dest, instantiate(unwrap_constmul(bc, Val(true)))) :
           copyto!(dest, instantiate(unwrap_constmul(bc, Val(false))))
end
false_only_copyto!(dest, bc) = copyto!(dest, instantiate(unwrap_constmul(bc, Val(false))))

lazy(c, a, b) = broadcasted(+, const_mul_broadcasted(LazyMulBy(typemax(Int)), c, a), b)
scalar(c, a, b) = broadcasted(+, const_mul_broadcasted(ScalarConstMul(), c, a), b)

f_C(out, c, a, b) = materialize!(out, scalar(c, a, b))
f_barrier(out, c, a, b) = materialize!(out, lazy(c, a, b))
f_inlined(out, c, a, b) = inlined_copyto!(out, lazy(c, a, b))
f_false(out, c, a, b) = false_only_copyto!(out, lazy(c, a, b))

const F = getfield(FiniteFieldArithmetic, Symbol(get(ARGS, 1, "GF2_16")))
const N = parse(Int, get(ARGS, 2, "4096"))

println("Julia ", VERSION, ", ", Sys.CPU_NAME, "; JULIA_LLVM_ARGS = '",
        get(ENV, "JULIA_LLVM_ARGS", ""), "'")
println(F, ", n = ", N, ", out .= c .* a .+ b below the threshold; median ns per element")
rng = Xoshiro(1)
a, b, out, c = rand(rng, F, N), rand(rng, F, N), Vector{F}(undef, N), rand(rng, F)
for (name, f) in (("C", f_C), ("barrier", f_barrier), ("inlined", f_inlined),
                  ("false only", f_false))
    t = @benchmark $f($out, $c, $a, $b) seconds = 1
    println("  ", rpad(name, 12), round(median(t).time / N; digits = 3), "  (allocs ", t.allocs, ")")
end

# Loop headers (blocks whose predecessors include themselves) and per-element `select`s of the
# broadcast index in the optimised IR.
bc = lazy(c, a, b)
for (name, f) in (("inlined", inlined_copyto!), ("false only", false_only_copyto!))
    ir = sprint(io -> code_llvm(io, f, (typeof(out), typeof(bc)); debuginfo = :none))
    loops = count(l -> (m = match(r"^(L[\w.]+):.*preds = (.*)$", l); m !== nothing &&
                        occursin("%" * m[1] * ",", m[2] * ",")), split(ir, '\n'))
    selects = count(l -> occursin(r"byteoffset[\w.]* = select", l), split(ir, '\n'))
    println("  IR of ", rpad(name, 12), loops, " loops, ", selects, " index selects")
end

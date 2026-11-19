# Runs the suite and prints the median time per element (or per step of the latency chain).
#   julia --project=bench bench/run.jl
# The bench environment's Manifest is not tracked; instantiate it once:
#   julia --project=bench -e 'using Pkg; Pkg.instantiate()'
include(joinpath(@__DIR__, "benchmarks.jl"))

tune!(SUITE)
results = run(SUITE; verbose = false, seconds = 2)
println("Julia ", VERSION, ", ", Sys.CPU_NAME, ", ", Sys.CPU_THREADS, " threads; elements per group: ",
        join(("$k: $v" for (k, v) in sort(collect(GROUP_SIZE))), ", "))
label(key) = join((k isa Tuple ? join(k, " / ") : k for k in key), " / ")
println(rpad("benchmark", 44), "median ns/element")
for (key, trial) in sort(collect(BenchmarkTools.leaves(results)); by = label ∘ first)
    println(rpad(label(key), 44), round(median(trial).time / GROUP_SIZE[first(key)]; digits = 3))
end

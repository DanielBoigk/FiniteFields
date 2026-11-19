# Profiles the multiplication strategies on the throughput workload of benchmarks.jl and prints
# the flat profile (Profile stdlib) and the allocations.
#   julia --project=bench bench/profile_nimber.jl
using Profile
include(joinpath(@__DIR__, "benchmarks.jl"))

rng = Xoshiro(1)
a, b = rand(rng, GF256, N), rand(rng, GF256, N)
c = similar(a)
for s in STRATEGIES
    f! = () -> broadcast!((x, y) -> mul(s, x, y), c, a, b)
    f!()                                    # compile
    println("== ", nameof(typeof(s)), ": allocated bytes per call = ", @allocated f!())
    Profile.clear()
    @profile for _ in 1:2000; f!(); end
    Profile.print(; format = :flat, sortedby = :count, mincount = 20, noisefloor = 2, C = false)
end

# Profiles the multiplication strategies on the throughput workload of benchmarks.jl and prints
# the allocations and the flat profile (Profile stdlib), for GF(256) and GF(2^128).
#   julia --project=bench bench/profile_nimber.jl
using Profile
include(joinpath(@__DIR__, "benchmarks.jl"))

rng = Xoshiro(1)
for (F, n, reps) in ((GF256, N, 2000), (GF2_128, N_LARGE, 500))
    a, b = rand(rng, F, n), rand(rng, F, n)
    c = similar(a)
    for s in STRATEGIES
        f! = () -> broadcast!((x, y) -> mul(s, x, y), c, a, b)
        f!()                                    # compile
        println("== ", F, " ", nameof(typeof(s)), ": allocated bytes per call = ", @allocated f!())
        Profile.clear()
        @profile for _ in 1:reps; f!(); end
        Profile.print(; format = :flat, sortedby = :count, mincount = 20, noisefloor = 2, C = false)
    end
end

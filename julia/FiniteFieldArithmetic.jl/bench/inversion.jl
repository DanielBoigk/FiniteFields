# J3: inversion strategies against each other and against the current default of `inv`, on
# nimber and polynomial-basis fields. Recorded in src/inversion.md.
# Run: julia --project=julia/FiniteFieldArithmetic.jl/bench julia/FiniteFieldArithmetic.jl/bench/inversion.jl
using FiniteFieldArithmetic, BenchmarkTools, Random
using FiniteFieldArithmetic: inverse, Fermat, ItohTsujii, Euclid, default_inv_strategy

const L = 4096

function inv_chain(s, x, y, n)          # latency: each inverse waits for the previous one
    for _ in 1:n
        x = inverse(s, x) + y
        iszero(x) && (x = y)
    end
    return x
end

ns(t) = round(t / L * 1e9; digits = 1)

const FIELDS = Any[("nimber GF(2^$(degree(F)))", F) for F in (GF256, GF2_16, GF2_32, GF2_64, GF2_128)]
append!(FIELDS, [("Conway GF(2^$n)", conway_field(n)) for n in (8, 16, 32, 33, 64, 65, 127)])
append!(FIELDS, [("AES GF(2^8)", GF256_AES), ("GHASH GF(2^128)", GF2_128_GHASH)])

println("Julia ", VERSION, ", ", Sys.CPU_NAME, "; ns per inverse, throughput / latency over $L elements")
println(rpad("field", 22), rpad("current default", 26), rpad("Fermat", 18), rpad("ItohTsujii", 18), "Euclid")
for (name, F) in FIELDS
    rng = Xoshiro(1)
    b = [x for x in rand(rng, F, 2L) if !iszero(x)][1:L]
    c = similar(b)
    cells = String[]
    for s in (default_inv_strategy(F), Fermat(), ItohTsujii(), Euclid())
        t = @belapsed $c .= inverse.(Ref($s), $b)
        l = @belapsed inv_chain($s, $(b[1]), $(b[2]), L)
        push!(cells, string(ns(t), " / ", ns(l)))
    end
    cells[1] = string(nameof(typeof(default_inv_strategy(F))), " ", cells[1])
    println(rpad(name, 22), rpad(cells[1], 26), rpad(cells[2], 18), rpad(cells[3], 18), cells[4])
end

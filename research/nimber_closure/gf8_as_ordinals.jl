# GF(8) inside Conway's On_2: its elements as ordinals. GF(8) ⊂ GF(64) = GF(4)(ω), ω^3 = 2
# (DiMuro Thm 3.1; d(ω) = 6, so ω ∉ GF(8)). An element a + bω + cω^2 (a, b, c ∈ GF(4), nimbers)
# is the ordinal [ω^2·c + ω·b + a] (DiMuro Lemma 1.1), ordered by (c, b, a) lexicographically.
# Run: julia --project=julia/FiniteFieldArithmetic.jl research/nimber_closure/gf8_as_ordinals.jl
using FiniteFieldArithmetic
mul3(x, y) = begin                       # (a0 + a1 ω + a2 ω^2)(b0 + …) mod ω^3 = 2
    r = [zero(GF4) for _ in 1:5]
    for i in 1:3, j in 1:3; r[i + j - 1] += x[i] * y[j]; end
    [r[1] + GF4(2) * r[4], r[2] + GF4(2) * r[5], r[3]]
end
pow3(x, e) = (r = [one(GF4), zero(GF4), zero(GF4)]; for _ in 1:e; r = mul3(r, x); end; r)
els = [[GF4(a), GF4(b), GF4(c)] for c in 0:3 for b in 0:3 for a in 0:3]   # ordinal order
f8 = filter(x -> pow3(x, 8) == x, els)
ordinal(x) = join(filter(!isempty, [Int(x[3].bits) == 0 ? "" : "ω²·$(Int(x[3].bits))",
                                    Int(x[2].bits) == 0 ? "" : "ω·$(Int(x[2].bits))",
                                    Int(x[1].bits) == 0 ? "" : string(Int(x[1].bits))]), " + ")
println("GF(8) in ordinal order (rank: ordinal):")
for (k, x) in enumerate(f8)
    println("  ", k - 1, ": ", isempty(ordinal(x)) ? "0" : ordinal(x))
end

# R1: compute Lenstra's α_u in On_2 for every odd prime u ≤ 127, from scratch, with On2.jl.
#
# For each prime u in increasing order (DiMuro, "On On_p", arXiv:1108.0962v3):
#   f = ord_u(2); χ_f = Σ_{r ∈ Q(f)} χ_r by Theorem 4.3 (u0 the smallest prime of h, r its power
#   in h, g = h/r: χ_h = χ_g if r | d(χ_g), else χ_g + χ_r), with the degrees d computed in towers;
#   α_u = [χ_f + m], m the excess, the first ordinal [χ_f + j], j = 0, 1, …, that is not a u-th
#   power in χ_u (Theorem 4.6; ordinals below χ_f are u-th powers, their degree not being
#   divisible by f). Then χ_u^u = α_u, d(χ_u) = u·d(α_u), and χ_{u^k}^u = χ_{u^(k-1)}.
# Lenstra's table (DiMuro Table 1, u ≤ 43) is recomputed, not used, and compared at the end.
#
# Run: julia --project=julia/FiniteFieldArithmetic.jl research/nimber_closure/compute_excesses.jl

include(joinpath(@__DIR__, "On2.jl"))
using .On2
using FiniteFieldArithmetic
using Printf

const LENSTRA = Dict(3 => 0, 5 => 0, 7 => 1, 11 => 1, 13 => 0, 17 => 0, 19 => 4, 23 => 1, 29 => 0,
                     31 => 1, 37 => 0, 41 => 1, 43 => 1)     # excess, DiMuro Table 1
const LENSTRA_Q = Dict(3 => [2], 5 => [4], 7 => [3], 11 => [5], 13 => [3, 4], 17 => [8],
                       19 => [9], 23 => [11], 29 => [4, 7], 31 => [5], 37 => [4, 9], 41 => [5],
                       43 => [7])

isprime(n) = n >= 2 && all(n % p != 0 for p in 2:isqrt(n))
ord2(u) = findfirst(k -> powermod(2, k, u) == 1, 1:u)
primepower(r) = (u = findfirst(p -> r % p == 0, 2:r) + 1; k = round(Int, log(u, r)); (u, k))

# What is known about the odd generators: r => (u, list of odd generators in the relation, scalar)
const REL = Dict{Int, Tuple{Int, Vector{Int}, Int}}()
const DEG = Dict{Int, Int}()        # d(χ_r) for prime powers r (2-powers: d(χ_{2^k}) = 2^k)

degchi(r) = ispow2(r) ? r : DEG[r]
twopow_value(r) = 1 << (r ÷ 2)      # χ_{2^k} = 2^(2^(k-1)) as a finite nimber

# Close a set of odd generators under their relations, in ordinal order.
function closure(rs)
    todo, seen = collect(rs), Set{Int}()
    while !isempty(todo)
        r = pop!(todo)
        r in seen && continue
        push!(seen, r)
        append!(todo, REL[r][2])
    end
    return sort!(collect(seen); by = r -> primepower(r))
end

# The tower over base B with the odd generators rs (closed), relations from REL.
function build_tower(B, rs)
    gens = Gen{B}[]
    for r in closure(rs)
        u, deps, c = REL[r]
        T = Tower(copy(gens))
        rel = element(T)
        for d in deps
            rel .+= gen_element(T, d)
        end
        rel[1] += B(c)
        push!(gens, Gen{B}(r, u, rel))
    end
    return Tower(gens)
end

# χ_h as (odd generators, finite scalar)
function chi_h(h)
    Q = Qset(h)
    odd = [r for r in Q if isodd(r)]
    c = sum((twopow_value(r) for r in Q if iseven(r)); init = 0)
    return odd, c
end

# ON2_BASE=GF2_16 runs everything over a larger base field (a check: results must not change)
const MINBASE = get(ENV, "ON2_BASE", "GF256") == "GF2_16" ? 65536 : 256
base_for(c) = max(c + 1, MINBASE) <= 256 ? GF256 : max(c + 1, MINBASE) <= 65536 ? GF2_16 : GF2_32

# degree of Σ χ_r over r in odd, plus the scalar c
function degree_of(odd, c)
    B = base_for(c)
    T = build_tower(B, odd)
    x = scalar_element(T, c)
    for r in odd
        x .+= gen_element(T, r)
    end
    return field_degree(T, x)
end

const QCACHE = Dict{Int, Vector{Int}}(1 => Int[])
function Qset(h)
    haskey(QCACHE, h) && return QCACHE[h]
    u0 = findfirst(p -> h % p == 0, 2:h) + 1
    r = u0
    while h % (r * u0) == 0
        r *= u0
    end
    g = h ÷ r
    qg = Qset(g)
    odd = [x for x in qg if isodd(x)]
    c = sum((twopow_value(x) for x in qg if iseven(x)); init = 0)
    dg = isempty(qg) ? 1 : degree_of(odd, c)
    Q = dg % r == 0 ? qg : sort!(vcat(qg, r))
    QCACHE[h] = Q
    return Q
end

results = []
println(" u   f(u)  Q(f)          excess  α_u                         d(α_u)    d(χ_u)   tower degree  time")
for u in 3:2:127
    isprime(u) || continue
    t0 = time()
    f = ord2(u)
    odd, c = chi_h(f)
    B = base_for(c + 64)                                  # room for c + j
    T = build_tower(B, odd)
    cache = Dict()
    χf = element(T)
    for r in odd
        χf .+= gen_element(T, r)
    end
    m = 0
    while true
        β = copy(χf)
        β[1] += B(c + m)                                  # [χ_f + m] = odd part ⊕ (c + m)
        is_upower(T, β, u, f, cache) || break
        m += 1
        m > 60 && error("excess search for u = $u did not stop")
    end
    α = copy(χf)
    α[1] += B(c + m)
    dα = field_degree(T, α, cache)
    REL[u] = (u, odd, c + m)
    DEG[u] = u * dα
    # prime powers u^k ≤ 128: χ_{u^k}^u = χ_{u^(k-1)}
    let r = u
        while r * u <= 128
            REL[r * u] = (u, [r], 0)
            DEG[r * u] = u * DEG[r]
            r *= u
        end
    end
    push!(results, (u, f, Qset(f), m, ordinal_string(T, α), dα, DEG[u], tower_degree(T)))
    @printf("%3d  %4d  %-12s  %5d   %-26s  %8d  %8d  %8d  %6.1fs\n", u, f, string(Qset(f)), m,
            ordinal_string(T, α), dα, DEG[u], tower_degree(T), time() - t0)
    flush(stdout)
end

println()
ok = true
for (u, f, Q, m, _, _, _, _) in results
    haskey(LENSTRA, u) || continue
    same = m == LENSTRA[u] && sort(Q) == sort(LENSTRA_Q[u])
    global ok &= same
    same || println("MISMATCH with Lenstra's table at u = $u: excess $m vs $(LENSTRA[u]), Q $Q vs $(LENSTRA_Q[u])")
end
println(ok ? "Lenstra's table (u ≤ 43) reproduced: excess and Q(f(u)) agree for all 13 primes." : "TABLE MISMATCH")

MINBASE == 256 && open(joinpath(@__DIR__, "data", "excesses.tsv"), "w") do io
    println(io, "# u\tf(u)\tQ(f(u))\texcess\talpha_u\td(alpha_u)\td(chi_u)")
    for (u, f, Q, m, α, dα, dχ, _) in results
        println(io, join((u, f, join(Q, ","), m, α, dα, dχ), '\t'))
    end
end

# The relations χ_r^u = (Σ odd generators) + scalar and the degrees d(χ_r), for rank_encoding.jl.
MINBASE == 256 && open(joinpath(@__DIR__, "data", "relations.jl"), "w") do io
    println(io, "# Generated by compute_excesses.jl: χ_r^u = Σ_{r' ∈ deps} χ_r' + scalar (scalar a finite")
    println(io, "# nimber), for the odd prime powers r ≤ 128 and the primes u ≤ 127; and d(χ_r) over GF(2).")
    println(io, "const RELATIONS = Dict{Int, Tuple{Int, Vector{Int}, Int}}(")
    for r in sort!(collect(keys(REL)))
        u, deps, c = REL[r]
        println(io, "    $r => ($u, $(repr(deps)), $c),")
    end
    println(io, ")")
    println(io, "const DEGREES = Dict{Int, Int}(", join(("$r => $(DEG[r])" for r in sort!(collect(keys(DEG)))), ", "), ")")
end

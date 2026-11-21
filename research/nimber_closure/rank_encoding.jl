# R1: the rank encoding of GF(2^n) ⊂ On_2 for n = 1 … 128 (ADR-0005).
#
# For each n: a tower L_n ⊇ GF(2^n) (base nimber field containing GF(2^(2-part of n)) and the
# relation scalars; odd generators χ_{u^k} for the odd prime powers u^k ∥ n, closed under their
# relations, from data/relations.jl); a spanning set of GF(2^n): products of a basis of each
# part (nimber bits 2^b for the 2-part; the Frobenius conjugates of θ_r = Tr_{L_n/GF(2^r)}(y) for
# an odd part r); its reduced echelon basis b_0 < … < b_{n-1} with respect to the ordinal order
# (highest digit first): rank(Σ ε_i b_i) = Σ ε_i 2^i, read off the pivot digits.
#
# Output per n (data/rank_encoding.jl), all in rank coordinates (bit i = coefficient of b_i):
#   g      the generator used internally: the basis element b_j of smallest j with degree n;
#   low    the minimal polynomial of g: x^n + Σ low_i x^i;
#   C      C[i+1] = rank(g^i), i = 0 … n-1 (the change of basis to the polynomial basis of g);
#   embed  for m | n, m < n: rank_n(b^(m)_i), i = 0 … m-1 (the inclusion GF(2^m) ⊂ GF(2^n)).
# The choice of g, y and of the spanning set does not affect any rank: they are canonical; only g
# (and so low, C) is an internal choice, documented as "the first basis element of degree n".
#
# Run: julia --project=julia/FiniteFieldArithmetic.jl research/nimber_closure/rank_encoding.jl [n1 n2 …]

include(joinpath(@__DIR__, "On2.jl"))
include(joinpath(@__DIR__, "data", "relations.jl"))
using .On2, FiniteFieldArithmetic, Random, Printf, Serialization
using FiniteFieldArithmetic: bits

const NIMBERS = (GF2, GF4, GF16, GF256, GF2_16, GF2_32, GF2_64, GF2_128)

primepower(r) = (u = findfirst(p -> r % p == 0, 2:r) + 1; (u, round(Int, log(u, r))))

function closure(rs)
    todo, seen = collect(rs), Set{Int}()
    while !isempty(todo)
        r = pop!(todo)
        r in seen && continue
        push!(seen, r)
        append!(todo, RELATIONS[r][2])
    end
    sort!(collect(seen); by = primepower)
end

function parts(n)
    a = 1 << trailing_zeros(n)                     # 2-part
    m, odd = n ÷ a, Int[]
    p = 3
    while m > 1
        if m % p == 0
            r = 1
            while m % p == 0
                m ÷= p; r *= p
            end
            push!(odd, r)
        end
        p += 2
    end
    a, odd
end

function base_for(n)
    a, odd = parts(n)
    need = max(a, 1)
    for r in closure(odd)
        c = RELATIONS[r][3]
        c > 0 && (need = max(need, 8 * sizeof(Int) - leading_zeros(c)))    # bits of the scalar
    end
    NIMBERS[findfirst(F -> degree(F) >= need, NIMBERS)]
end

function build_tower(B, odd)
    gens = Gen{B}[]
    for r in closure(odd)
        u, deps, c = RELATIONS[r]
        T = Tower(copy(gens))
        rel = element(T)
        for d in deps
            rel .+= gen_element(T, d)
        end
        rel[1] += B(c)
        push!(gens, Gen{B}(r, u, rel))
    end
    Tower(gens)
end

# --- Digits: an element as a bit string, highest digit = most significant in ordinal order ----

w(T) = degree(eltype(element(T)))
digit(T, x, k) = isodd(bits(x[k ÷ w(T) + 1]) >> (k % w(T)))      # digit k (0-based)
function topdigit(T, x)
    for idx in length(x):-1:1
        b = bits(x[idx])
        iszero(b) && continue
        return (idx - 1) * w(T) + (8 * sizeof(b) - 1 - leading_zeros(b))
    end
    return -1
end

# Reduced echelon basis of linearly independent elements, sorted by pivot ascending.
function rref(T, vs)
    rows = Tuple{Int, Vector}[]                   # (pivot, element)
    for v in vs
        x = copy(v)
        for (p, r) in sort(rows; by = first, rev = true)
            digit(T, x, p) && (x = x .+ r)
        end
        p = topdigit(T, x)
        p < 0 && error("dependent vectors")
        # clear the new pivot from the existing rows
        rows = [(q, digit(T, r, p) ? r .+ x : r) for (q, r) in rows]
        push!(rows, (p, x))
    end
    sort!(rows; by = first)
end

"""Rank coordinates of x in the echelon basis; checks that x lies in the span."""
function rankof(T, basis, x)
    k = zero(UInt128)
    y = copy(x)
    for (i, (p, b)) in enumerate(basis)
        if digit(T, x, p)
            k |= UInt128(1) << (i - 1)
            y = y .+ b
        end
    end
    all(iszero, y) || error("element not in the subfield")
    k
end

element_of(T, basis, k::UInt128) = (x = element(T); for i in 0:(length(basis) - 1)
    isodd(k >> i) && (x = x .+ basis[i + 1][2]); end; x)

# --- One n ------------------------------------------------------------------------------------

function compute(n; verbose = true)
    t0 = time()
    a, odd = parts(n)
    B = base_for(n)
    T = build_tower(B, odd)
    cache = Dict()
    rng = Xoshiro(n)
    # part bases
    pbases = Vector{Vector{Vector{B}}}()
    push!(pbases, [scalar_element(T, one(UInt128) << b) for b in 0:(a - 1)])   # nimber bits of GF(2^a)
    for r in odd
        F1 = frob1(T)
        while true
            y = [rand(rng, B) for _ in 1:T.S]
            θ = trace_to(T, y, r, cache)
            conj = [θ]
            for _ in 2:r
                push!(conj, apply(T, F1, conj[end]))
            end
            apply(T, F1, conj[end]) == θ || error("θ is not in GF(2^$r)")
            try
                rref(T, conj)
                push!(pbases, conj)
                break
            catch
                verbose && println("  n = $n: θ_$r not normal, another y")
            end
        end
    end
    # spanning set: products, one factor from each part
    span = [scalar_element(T, 1)]
    for pb in pbases
        span = [length(pb) == 1 ? s : (pb === pbases[1] ? s .* Ref(p[1]) : mul(T, s, p)) for s in span for p in pb]
    end
    length(span) == n || error("spanning set has $(length(span)) elements, not $n")
    basis = rref(T, span)
    rankof(T, basis, scalar_element(T, 1)) == 1 || error("b_0 ≠ 1")
    # internal generator: smallest rank with degree n
    # internal generator g: the basis element b_j with the smallest j whose degree is n (then, if
    # none, ranks 2^(n-1) + i); an implementation choice, no rank depends on it
    g, gk = nothing, UInt128(0)
    n == 1 && ((g, gk) = (scalar_element(T, 1), UInt128(1)))
    cands = n == 1 ? UInt128[] : vcat([UInt128(1) << j for j in 1:(n - 1)],
                                      [(UInt128(1) << (n - 1)) + i for i in 1:1000])
    for k in cands
        x = element_of(T, basis, k)
        if field_degree(T, x, cache) % n == 0
            g, gk = x, k
            break
        end
    end
    g === nothing && error("no generator found")
    # multiplication by g on the rank basis, then powers and minimal polynomial (in coordinates)
    Mg = [rankof(T, basis, mul(T, g, b)) for (_, b) in basis]
    matvec(M, v) = (r = zero(UInt128); for i in 0:(n - 1); isodd(v >> i) && (r ⊻= M[i + 1]); end; r)
    C = UInt128[1]
    for _ in 2:n
        push!(C, matvec(Mg, C[end]))
    end
    gn = matvec(Mg, C[end])                         # g^n
    low = solve(C, gn, n)                           # g^n = Σ low_i g^i
    verbose && @printf("n = %3d  tower %-24s base GF(2^%d)  degree %7d  g = rank %d  %.1fs\n", n,
                       join(("χ$(g.r)" for g in T.gens), ","), degree(B), tower_degree(T), gk, time() - t0)
    (; n, T, basis, gk, low, C, Mg)
end

# solve Σ c_i C[i] = v over GF(2) (C a basis), return c as bits
function solve(C, v, n)
    rows = [(C[i], UInt128(1) << (i - 1)) for i in 1:n]
    piv = Tuple{UInt128, UInt128}[]
    for (r, c) in rows
        for (pr, pc) in piv
            if !iszero(r & (UInt128(1) << (127 - leading_zeros(pr))))
                r ⊻= pr; c ⊻= pc
            end
        end
        iszero(r) && error("powers of g are dependent: degree < n")
        push!(piv, (r, c))
        sort!(piv; by = x -> leading_zeros(x[1]))
    end
    c = zero(UInt128)
    for (pr, pc) in piv
        if !iszero(v & (UInt128(1) << (127 - leading_zeros(pr))))
            v ⊻= pr; c ⊻= pc
        end
    end
    iszero(v) || error("not in the span")
    c
end

# Map an element of the tower Tm into the tower Tn ⊇ Tm (generators and base included).
function lift(Tm, Tn, x)
    Bn = eltype(element(Tn))
    y = element(Tn)
    jmap = [findfirst(g -> g.r == gm.r, Tn.gens) for gm in Tm.gens]
    any(isnothing, jmap) && error("tower not contained")
    for idx in 1:length(x)
        iszero(x[idx]) && continue
        e = idx - 1
        f = 0
        for (j, jn) in enumerate(jmap)
            ej = (e ÷ Tm.s[j]) % Tm.u[j]
            f += ej * Tn.s[jn]
        end
        y[f + 1] += Bn(bits(x[idx]))
    end
    y
end

ns = isempty(ARGS) ? collect(1:128) : parse.(Int, ARGS)
outdir = joinpath(@__DIR__, "data", "rank")
mkpath(outdir)
results = Dict{Int, Any}()
for n in sort(ns)
    R = compute(n)
    # embeddings from every proper divisor m (recomputed or loaded)
    emb = Dict{Int, Vector{UInt128}}()
    for m in 1:(n - 1)
        n % m == 0 || continue
        Rm = get(results, m, nothing)
        if Rm === nothing
            f = joinpath(outdir, "basis_$m.jls")
            Rm = isfile(f) ? deserialize(f) : compute(m; verbose = false)
        end
        emb[m] = [rankof(R.T, R.basis, lift(Rm.T, R.T, b)) for (_, b) in Rm.basis]
    end
    open(joinpath(outdir, "n_$n.txt"), "w") do io
        println(io, "n\t$n")
        println(io, "g\t$(R.gk)")
        println(io, "low\t$(R.low)")
        println(io, "C\t", join(R.C, ","))
        for m in sort!(collect(keys(emb)))
            println(io, "embed $m\t", join(emb[m], ","))
        end
    end
    serialize(joinpath(outdir, "basis_$n.jls"), (; T = R.T, basis = R.basis))
    results[n] = (; T = R.T, basis = R.basis)
    flush(stdout)
end

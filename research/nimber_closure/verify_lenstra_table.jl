# Checks Table 1 of DiMuro, "On On_p" (arXiv:1108.0962v3), the elements α_u of Conway's On_2 for
# odd primes u ≤ 43 (first printed in Lenstra, "On the algebraic closure of two", 1977), with
# the package's nimber arithmetic and a generic radical extension x^u - a.
#
# What is checked, per u (DiMuro Thm 3.1(4), Thm 4.6): α_u = [χ_f + m] (f = ord_u 2, m the
# "excess") is not a u-th power in the field χ_u, and every ordinal [χ_f + j], j < m, is. Ordinals
# below χ_f are u-th powers automatically (their degree d is not divisible by f, so u ∤ 2^d - 1).
#
# Test for β in a finite field F of degree E with u | 2^E - 1 and u ∤ E / d(β) (true here: E is a
# product of 2-powers and primes < u): β is a u-th power in χ_u iff β^((2^E - 1)/u) = 1
# (DiMuro Lemma 3.6; extending F by degrees prime to u does not change the u-part of 2^E - 1).
#
# Fields used (ordinals as nim-polynomials, DiMuro Lemma 1.1): χ_2 = 2, χ_4 = 4, χ_8 = 16 (finite
# nimbers); χ_3 = ω with ω^3 = 2; χ_9 = ω^3 (ordinal) with χ_9^3 = ω; χ_5 = ω^ω with χ_5^5 = 4;
# χ_7 = ω^(ω^2) with χ_7^7 = ω + 1; χ_11 = ω^(ω^3) with χ_11^11 = χ_5 + 1 (Thm 3.1(4), (5)).
#
# Run: julia --project=julia/FiniteFieldArithmetic.jl research/nimber_closure/verify_lenstra_table.jl

using FiniteFieldArithmetic

# --- A radical extension B[x]/(x^n - a), elements as coefficient vectors -----------------------

struct Ext{B}
    n::Int
    a::B
    name::String
end
struct El{B}
    K::Ext{B}
    c::Vector{B}
end
Base.zero(K::Ext{B}) where {B} = El(K, [zero_of(K.a) for _ in 1:K.n])
Base.one(K::Ext{B}) where {B} = (z = zero(K); z.c[1] = one_of(K.a); z)
gen(K::Ext) = (z = zero(K); z.c[2] = one_of(K.a); z)
lift(K::Ext{B}, b::B) where {B} = (z = zero(K); z.c[1] = b; z)
zero_of(b) = zero(b)
one_of(b) = one(b)
zero_of(e::El) = zero(e.K)
one_of(e::El) = one(e.K)
Base.zero(e::El) = zero(e.K)
Base.one(e::El) = one(e.K)
Base.:+(x::El, y::El) = El(x.K, x.c .+ y.c)
Base.:(==)(x::El, y::El) = x.c == y.c
Base.isone(x::El) = x == one(x.K)
function Base.:*(x::El{B}, y::El{B}) where {B}
    K, n = x.K, x.K.n
    r = [zero_of(K.a) for _ in 1:(2n - 1)]
    for i in 1:n, j in 1:n
        r[i + j - 1] = r[i + j - 1] + x.c[i] * y.c[j]
    end
    for k in (2n - 1):-1:(n + 1)              # x^(n + t) = a · x^t
        r[k - n] = r[k - n] + K.a * r[k]
    end
    return El(K, r[1:n])
end
function pow(x, e::Integer)
    r, b = one_of(x), x
    while e > 0
        isodd(e) && (r = r * b)
        b = b * b
        e >>= 1
    end
    return r
end

# u-th power test in a field of degree E over GF(2)
is_upower(β, u, E) = isone(pow(β, (big(2)^E - 1) ÷ u))

# --- The fields ------------------------------------------------------------------------------

const N = GF256                                   # finite nimbers up to 255 (contains 2, 4, 16)
n(x) = N(x)
K3 = Ext(3, n(2), "ω")                            # F(ω), ω^3 = 2
ω = gen(K3)
K9 = Ext(3, ω, "χ9")                              # F(ω)(χ_9), χ_9^3 = ω
K5 = Ext(5, n(4), "χ5")                           # F(χ_5), χ_5^5 = 4
χ5 = gen(K5)
K7 = Ext(7, ω + one(K3), "χ7")                    # F(ω)(χ_7), χ_7^7 = ω + 1
K11 = Ext(11, χ5 + one(K5), "χ11")                # F(χ_5)(χ_11), χ_11^11 = χ_5 + 1

# β(j) = [χ_f + j] as a field element, and the degree E of the field it is computed in.
# (For a group χ ≥ ω and finite j, the ordinal sum [χ + j] is the nim-sum χ + j; for χ_8 = 16 and
# j < 16 too; for χ_12 = [ω + 4]: [ω + 4 + j] = ω + (4 + j) with 4 + j an ordinary integer.)
const CASES = [
    # u, f, excess, β(j), E
    (3, 2, 0, j -> n(2 + j), 8),
    (5, 4, 0, j -> n(4 + j), 8),
    (7, 3, 1, j -> ω + lift(K3, n(j)), 24),
    (11, 10, 1, j -> χ5 + lift(K5, n(j)), 40),
    (13, 12, 0, j -> ω + lift(K3, n(4 + j)), 24),
    (17, 8, 0, j -> n(16 + j), 8),
    (19, 18, 4, j -> gen(K9) + lift(K9, lift(K3, n(j))), 72),
    (23, 11, 1, j -> gen(K11) + lift(K11, lift(K5, n(j))), 440),
    (29, 28, 0, j -> gen(K7) + lift(K7, lift(K3, n(4 + j))), 168),
    (31, 5, 1, j -> χ5 + lift(K5, n(j)), 40),
    (37, 36, 0, j -> gen(K9) + lift(K9, lift(K3, n(4 + j))), 72),
    (41, 20, 1, j -> χ5 + lift(K5, n(j)), 40),
    (43, 14, 1, j -> gen(K7) + lift(K7, lift(K3, n(j))), 168),
]

println("u   f   excess  [χ_f + j] is a u-th power for j = 0..excess   (expected: yes … yes, no)")
allok = true
for (u, f, m, β, E) in CASES
    @assert (big(2)^E - 1) % u == 0
    @assert f == findfirst(k -> (big(2)^k - 1) % u == 0, 1:u)         # f = ord_u 2
    res = [is_upower(β(j), u, E) for j in 0:m]
    ok = all(res[1:m]) && !res[m + 1]
    global allok &= ok
    println(rpad(u, 4), rpad(f, 4), rpad(m, 8), join((r ? "yes" : "no") for r in res), "   ", ok ? "ok" : "MISMATCH")
end
# the defining relations are not u-th powers either (the extensions are fields: Abel–Capelli)
for (x, u, E, s) in ((n(2), 3, 8, "2 (ω^3)"), (ω, 3, 24, "ω (χ_9^3)"), (n(4), 5, 8, "4 (χ_5^5)"),
                     (ω + one(K3), 7, 24, "ω+1 (χ_7^7)"), (χ5 + one(K5), 11, 40, "χ_5+1 (χ_11^11)"))
    ok = !is_upower(x, u, E)
    global allok &= ok
    println("relation ", rpad(s, 18), " not a $u-th power: ", ok)
end
println(allok ? "ALL CHECKS PASS" : "SOME CHECKS FAIL")

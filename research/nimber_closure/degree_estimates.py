# Estimates of d(χ_r), the degree over GF(2) of Conway's canonical generators χ_r of On_2, for
# prime powers r with primes ≤ 127; these are the field degrees in which α_u (and χ_u) must be
# computed. Structure from DiMuro, "On On_p", Thm 3.1, Thm 4.3, Thm 4.6:
#   d(χ_{2^k}) = 2^k;  χ_u^u = α_u, so d(χ_u) = u · d(α_u);  χ_{u^{k+1}}^u = χ_{u^k}, so
#   d(χ_{u^{k+1}}) = u · d(χ_{u^k});  α_u = [χ_f + m], f = ord_u(2), m the excess;
#   χ_h by Thm 4.3: u0 smallest prime of h, r its power in h, g = h / r:
#   χ_h = χ_g if r | d(χ_g), else χ_g + χ_r.
# Assumptions (ESTIMATE, not proved): d(a + b) = lcm(d(a), d(b)) for the sums occurring here, and
# for u > 43 (excess unknown) d(m) ≤ 4, as for every u ≤ 43 (Table 1: excess ≤ 4).
from math import lcm
from sympy import factorint, primerange, n_order

KNOWN_EXCESS = {3:0,5:0,7:1,11:1,13:0,17:0,19:4,23:1,29:0,31:1,37:0,41:1,43:1}
def dnim(m):            # degree of the finite nimber m over GF(2): 2^level
    if m <= 1: return 1
    k = 0
    while m >= 2 ** (2 ** k): k += 1
    return 2 ** k

dchi = {}               # d(χ_r) for prime powers r
def d_chi_r(u, k):
    r = u ** k
    if r in dchi: return dchi[r]
    if u == 2: v = r
    elif k > 1: v = u * d_chi_r(u, k - 1)
    else:
        f = n_order(2, u)
        m = KNOWN_EXCESS.get(u)
        v = u * lcm(d_chi_h(f), dnim(m) if m is not None else 4)
    dchi[r] = v
    return v

def Q(h):               # the prime powers of χ_h = Σ_{r∈Q} χ_r (Thm 4.3)
    if h == 1: return []
    fs = sorted(factorint(h).items())
    u0, k0 = fs[0]
    r, g = u0 ** k0, h // u0 ** k0
    qg = Q(g)
    dg = lcm(1, *[dchi_of(x) for x in qg]) if qg else 1
    return qg if dg % r == 0 else qg + [r]
def dchi_of(r):
    (u, k), = factorint(r).items()
    return d_chi_r(u, k)
def d_chi_h(h):
    return lcm(1, *[dchi_of(r) for r in Q(h)])

print("u    f(u)  Q(f)        d(χ_u) (estimate)   digits of d")
for u in primerange(3, 128):
    f = n_order(2, u)
    d = d_chi_r(u, 1)
    flag = "" if u in KNOWN_EXCESS else "  (excess unknown)"
    print(f"{u:<4} {f:<5} {str(sorted(Q(f))):<11} {d:<20} {len(str(d))}{flag}")

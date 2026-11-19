# Cost of the parts of a CLMul product (src/nimber_fast.md, "Breakdown of CLMul"): the change of
# basis `to_poly`, the product in the polynomial basis `poly_mulmod`, and the whole product on bit
# patterns, per element over 4096 random values.
#   julia --project=bench bench/clmul_parts.jl
using BenchmarkTools, FiniteFieldArithmetic, Random
const M = FiniteFieldArithmetic

M.has_pclmul() || @warn "no PCLMULQDQ: poly_mulmod uses the software carry-less product"
for (T, V) in ((UInt64, Val(6)), (UInt128, Val(7)))
    rng = Xoshiro(1)
    a, b = rand(rng, T, 4096), rand(rng, T, 4096)
    c = similar(a)
    full = (x, y) -> M.nim_mul_bits(M.CLMul(), x, y, V)
    t_to = @belapsed broadcast!(M.to_poly, $c, $a)
    t_mm = @belapsed broadcast!(M.poly_mulmod, $c, $a, $b)
    t_full = @belapsed broadcast!($full, $c, $a, $b)
    t_tab = @belapsed broadcast!((x, y) -> M.nim_mul_bits(M.MulTable(), x, y, $V), $c, $a, $b)
    ns(t) = round(t / 4096 * 1e9; digits = 2)
    println(T, ": to_poly ", ns(t_to), " ns, poly_mulmod ", ns(t_mm), " ns, CLMul product ",
            ns(t_full), " ns, MulTable product ", ns(t_tab), " ns")
end

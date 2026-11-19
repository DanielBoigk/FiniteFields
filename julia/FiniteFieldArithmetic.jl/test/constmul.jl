# Multiplication of an array by a constant in broadcasts (backlog J8c): `c .* a` and `a .* c`
# go through `default_constmul_strategy`, which picks one of
# - `ScalarConstMul()` (C): the scalar `*` element by element (the behaviour before J8c);
# - `EagerMulBy()` (A): `MulBy(c).(a)` at once, which materialises inside fused expressions;
# - `LazyMulBy(t)` (B): a lazy broadcast of a callable that holds `c` and builds the `MulBy`
#   tables only for arguments of at least `t` elements, so fusion is kept.
#
# Oracle: `mul(Recursive(), …)`, independent of all three. The strategy is switched by redefining
# `default_constmul_strategy` for one field (the documented way), and every check then runs in
# the latest world through `invokelatest`.

using FiniteFieldArithmetic
using FiniteFieldArithmetic: bits, mul, Recursive, ConstMulStrategy, ScalarConstMul, EagerMulBy,
    LazyMulBy, ConstMul, const_mul_broadcasted, default_constmul_strategy, mulby_threshold
using Random, Test

rng = Xoshiro(20261122)
const CM_LEVELS = (GF2, GF4, GF16, GF256, GF2_16, GF2_32, GF2_64, GF2_128)
cm_refmul(a, b) = mul(Recursive(), a, b)

# Run `f()` with `default_constmul_strategy(F) == s`, then remove the redefinition again.
function with_constmul(f, ::Type{F}, s::ConstMulStrategy) where {F}
    @eval FiniteFieldArithmetic.default_constmul_strategy(::Type{$F}) = $s
    try
        return Base.invokelatest(f)
    finally
        Base.delete_method(which(FiniteFieldArithmetic.default_constmul_strategy, (Type{F},)))
    end
end

# The expressions under test, as functions so that allocations are measured on compiled code.
cm_left(c, a) = c .* a
cm_right(a, c) = a .* c
cm_fused(out, c, a, b) = (out .= c .* a .+ b; out)
cm_into(out, c, a) = (out .= c .* a; out)
cm_lazyarg(c, a, b) = c .* (a .+ b)
cm_nested(c, d, a) = d .* (c .* a)
cm_star(c, a) = c * a

function check_constmul(F, s, n)
    a, b = rand(rng, F, n), rand(rng, F, n)
    n >= 2 && (a[1] = zero(F); a[2] = one(F))
    c, d = rand(rng, F), rand(rng, F)
    ref = cm_refmul.(Ref(c), a)
    @test cm_left(c, a) == ref
    @test cm_left(c, a) isa Vector{F}
    @test cm_right(a, c) == ref
    @test cm_fused(similar(a), c, a, b) == ref .+ b
    @test cm_into(similar(a), c, a) == ref
    @test cm_lazyarg(c, a, b) == cm_refmul.(Ref(c), a .+ b)
    @test cm_nested(c, d, a) == cm_refmul.(Ref(d), ref)
    @test cm_star(c, a) == ref
    M = reshape(rand(rng, F, 3n), 3, n)
    @test cm_left(c, M) == cm_refmul.(Ref(c), M)              # any array shape
    @test cm_left(zero(F), a) == zeros(F, n)
    @test cm_left(one(F), a) == a
    return nothing
end

@testset "default is the scalar product" begin
    @test ScalarConstMul() isa ConstMulStrategy
    @test EagerMulBy() isa ConstMulStrategy
    @test LazyMulBy(10) isa ConstMulStrategy
    for F in CM_LEVELS
        @test default_constmul_strategy(F) === ScalarConstMul()
        @test mulby_threshold(F) isa Int && mulby_threshold(F) >= 0
        @test LazyMulBy(F) === LazyMulBy(mulby_threshold(F))
    end
end

@testset "c .* a with $(nameof(typeof(s))) $(s isa LazyMulBy ? s.threshold : ""), $F" for F in CM_LEVELS,
        s in (ScalarConstMul(), EagerMulBy(), LazyMulBy(0), LazyMulBy(typemax(Int)), LazyMulBy(50))
    with_constmul(F, s) do
        for n in (0, 1, 7, 49, 50, 300)
            check_constmul(F, s, n)
        end
    end
end

@testset "mixed and non-array operands are left alone, $F" for F in (GF16, GF2_64)
    with_constmul(F, LazyMulBy(0)) do
        a = rand(rng, F, 20)
        @test GF256(5) .* a == GF256(5) .* GF256.(a)             # other field: scalar product
        @test (rand(rng, F, 4) .* rand(rng, F)) isa Vector{F}
        x, y = rand(rng, F), rand(rng, F)
        @test x .* y == cm_refmul(x, y)                          # scalars: no arrays
        @test a .* a == cm_refmul.(a, a)                         # two arrays: not a constant
    end
end

@testset "lazy strategy: tables only above the threshold, $F" for F in CM_LEVELS
    c, a = rand(rng, F), rand(rng, F, 64)
    bc = const_mul_broadcasted(LazyMulBy(64), c, a)
    @test bc isa Base.Broadcast.Broadcasted
    m = bc.f
    @test m isa ConstMul{F}
    @test m.tables                                            # 64 ≥ 64
    @test !const_mul_broadcasted(LazyMulBy(65), c, a).f.tables
    @test typeof(const_mul_broadcasted(LazyMulBy(65), c, a)) == typeof(bc)   # type-stable
    lazy = Base.Broadcast.broadcasted(+, a, a)                # a lazy argument of length 64
    @test const_mul_broadcasted(LazyMulBy(64), c, lazy).f.tables
    @test !const_mul_broadcasted(LazyMulBy(65), c, lazy).f.tables
    for x in rand(rng, F, 50)
        @test m(x) == cm_refmul(c, x)
        @test const_mul_broadcasted(LazyMulBy(65), c, a).f(x) == cm_refmul(c, x)
    end
end

@testset "eager strategy materialises at once, $F" for F in (GF256, GF2_32)
    c, a = rand(rng, F), rand(rng, F, 10)
    r = const_mul_broadcasted(EagerMulBy(), c, a)
    @test r isa Vector{F} && r == cm_refmul.(Ref(c), a)
    @test Base.Broadcast.materialize(r) === r
end

# Allocations: the scalar strategy and the lazy one below the threshold allocate nothing in fused
# expressions; above the threshold the lazy one allocates its tables only (and up to GF(256)
# nothing, its table being a column of the product table).
@testset "allocations, $F" for F in CM_LEVELS
    a, b, out = rand(rng, F, 1000), rand(rng, F, 1000), Vector{F}(undef, 1000)
    c = rand(rng, F)
    tabbytes = degree(F) <= 8 ? 0 : 256 * sizeof(FiniteFieldArithmetic.storagetype(F))^2
    for (s, limit) in ((ScalarConstMul(), 0), (LazyMulBy(typemax(Int)), 0),
                       (LazyMulBy(0), tabbytes == 0 ? 0 : tabbytes + 256))
        with_constmul(F, s) do
            cm_fused(out, c, a, b); cm_into(out, c, a)           # compile
            @test (@allocated cm_fused(out, c, a, b)) <= limit
            @test (@allocated cm_into(out, c, a)) <= limit
        end
    end
end

# Multiplication of an array by a constant in broadcasts (backlog J8c, J8d): `c .* a` and
# `a .* c` go through `default_constmul_strategy`, which picks one of
# - `ScalarConstMul()`: the scalar `*` element by element (Base's behaviour; default up to GF(256));
# - `LazyMulBy(t)`: a lazy broadcast of `MulBy(c)` for arguments of at least `t` elements, and
#   Base's broadcast below, so fusion is kept (default from GF(2^16), `t = mulby_threshold(F)`).
#
# Oracle: `mul(Recursive(), …)`, independent of both. The strategy is switched by redefining
# `default_constmul_strategy` for one field (the documented way), and every check then runs in
# the latest world through `invokelatest`.

using FiniteFieldArithmetic
using FiniteFieldArithmetic: bits, mul, Recursive, ConstMulStrategy, ScalarConstMul, LazyMulBy,
    const_mul_broadcasted, default_constmul_strategy, mulby_threshold
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

@testset "defaults: scalar up to GF(256), lazy tables above" begin
    @test ScalarConstMul() isa ConstMulStrategy
    @test LazyMulBy(10) isa ConstMulStrategy
    for F in CM_LEVELS
        @test mulby_threshold(F) isa Int && mulby_threshold(F) >= 0
        @test LazyMulBy(F) === LazyMulBy(mulby_threshold(F))
        @test default_constmul_strategy(F) ===
              (degree(F) <= 8 ? ScalarConstMul() : LazyMulBy(F))
    end
    # the default reaches the dot syntax: long arrays use tables, short ones Base's broadcast
    c = GF2_32(0x53)
    long, short = rand(rng, GF2_32, mulby_threshold(GF2_32)), rand(rng, GF2_32, 3)
    @test Base.Broadcast.broadcasted(*, c, long).f isa MulBy
    @test Base.Broadcast.broadcasted(*, c, short).f === (*)
    @test c .* long == cm_refmul.(Ref(c), long)
end

@testset "c .* a with $(nameof(typeof(s))) $(s isa LazyMulBy ? s.threshold : ""), $F" for F in CM_LEVELS,
        s in (ScalarConstMul(), LazyMulBy(0), LazyMulBy(typemax(Int)), LazyMulBy(50))
    with_constmul(F, s) do
        for n in (0, 1, 7, 49, 50, 300)
            check_constmul(F, s, n)
        end
    end
end

@testset "mixed and non-array operands are left alone, $F" for F in (GF16, GF2_64)
    with_constmul(F, LazyMulBy(0)) do
        a = rand(rng, F, 20)
        @test GF2(1) .* a == a                                   # other field: scalar product
        @test GF4(2) .* a == cm_refmul.(Ref(F(2)), a)
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
    @test m isa MulBy{F}                                      # 64 ≥ 64: tables
    below = const_mul_broadcasted(LazyMulBy(65), c, a)
    @test below.f === (*) && below.args === (c, a)              # Base's own broadcast
    @test const_mul_broadcasted(LazyMulBy(65), a, c, (a, c)).args === (a, c)  # order kept
    # The result is one of two Broadcasted types, a small Union that the compiler splits.
    rt = Base.return_types(const_mul_broadcasted, Tuple{LazyMulBy, F, Vector{F}})
    @test length(rt) == 1 && rt[1] == Union{typeof(bc), typeof(below)}
    lazy = Base.Broadcast.broadcasted(+, a, a)                # a lazy argument of length 64
    @test const_mul_broadcasted(LazyMulBy(64), c, lazy).f isa MulBy{F}
    @test const_mul_broadcasted(LazyMulBy(65), c, lazy).f === (*)
    for x in rand(rng, F, 50)
        @test m(x) == cm_refmul(c, x)
    end
    @test m(GF2(1)) == c                                      # other operand types: scalar `*`
end

# Allocations: the scalar strategy and the lazy one below the threshold (a Union result, split by
# the compiler) allocate nothing in fused
# expressions; above the threshold the lazy one allocates its tables only (256 B at GF(256) …
# 64 KiB at GF(2^128)).
@testset "allocations, $F" for F in CM_LEVELS
    a, b, out = rand(rng, F, 1000), rand(rng, F, 1000), Vector{F}(undef, 1000)
    c = rand(rng, F)
    tabbytes = 256 * sizeof(FiniteFieldArithmetic.storagetype(F))^2
    for (s, limit) in ((ScalarConstMul(), 0), (LazyMulBy(typemax(Int)), 0),
                       (LazyMulBy(0), tabbytes + 256))
        with_constmul(F, s) do
            cm_fused(out, c, a, b); cm_into(out, c, a)           # compile
            @test (@allocated cm_fused(out, c, a, b)) <= limit
            @test (@allocated cm_into(out, c, a)) <= limit
        end
    end
end

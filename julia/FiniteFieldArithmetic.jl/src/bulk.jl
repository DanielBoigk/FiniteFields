# Bulk and array operations (backlog J8): uniform sampling of arrays, multiplication by a fixed
# constant, matrix products, and bitsliced vectors. Theory and measurements: bulk.md beside this
# file. Broadcasting needs no code of its own: the element types are `isbits` numbers, so
# `a .* b`, `c .* a`, `inv.(a)` and fused expressions compile to plain loops over the scalar
# arithmetic (tested in test/bulk.jl).

# --- Random sampling --------------------------------------------------------------------------

# Arrays are filled with random bit patterns in one call on the storage type (which the RNG
# does in bulk, with SIMD for Xoshiro), then masked when the field does not fill its storage.
# Uniform: every bit pattern of the storage type is equally likely, and the mask keeps the low
# n bits. The elements differ from those of per-element sampling with the same seed.
function Random.rand!(rng::Random.AbstractRNG, A::Array{F},
                      ::Random.SamplerType{F}) where {F <: BinaryFieldElement}
    T = storagetype(F)
    GC.@preserve A Random.rand!(rng, unsafe_wrap(Array, Ptr{T}(pointer(A)), size(A)))
    if !fills_storage(F)
        m = max_bits(F)
        @inbounds for i in eachindex(A)
            A[i] = unsafe_frombits(F, bits(A[i]) & m)
        end
    end
    return A
end

# --- Multiplication by a constant -------------------------------------------------------------

"""
    MulBy(c::Nimber)

Multiplication by the constant `c`, prepared for many operands: `m = MulBy(c)`, then `m(x)`
equals `c * x`, and `m.(a)` multiplies an array. ``x \\mapsto c \\otimes x`` is linear over
``\\mathbb{F}_2``, so it is stored as byte tables: `m(x)` is one lookup per byte of `x` and XORs
(``n/8`` lookups in ``\\mathbb{F}_{2^n}``, one for ``n \\le 8``). Setting up costs ``n``
products and ``32 n`` XORs and ``32 n`` storage words (64 KiB at ``\\mathbb{F}_{2^{128}}``), so
it pays for long arrays only; the measurements are in `src/bulk.md`.
"""
struct MulBy{N <: Nimber, T <: Unsigned} <: Function
    c::N
    tab::Matrix{T}
end

function MulBy(c::N) where {K, T, N <: Nimber{K, T}}
    tab = Matrix{T}(undef, 256, sizeof(T))
    return MulBy{N, T}(c, mulby_tables!(tab, c))
end

# The byte tables of x ↦ c ⊗ x, written into `tab` (256 × sizeof(T)).
function mulby_tables!(tab::Matrix{T}, c::Nimber{K, T}) where {K, T}
    if K <= 3                                     # GF(256) and below: a column of the table
        @inbounds for b in 1:256
            tab[b, 1] = NIM256_MUL[b, Int(bits(c)) + 1] % T
        end
        return tab
    end
    s = default_strategy(typeof(c))
    img = ntuple(k -> nim_mul_bits(s, bits(c), one(T) << (k - 1), Val(K)), Val(8 * sizeof(T)))
    return byte_tables!(tab, img)
end
byte_tables!(tab::Matrix{T}, img::NTuple{M, T}) where {M, T <: Unsigned} =
    byte_tables!(tab, StaticImages(img))

# A read-only vector view of a tuple, so byte_tables! can take the images without allocating.
struct StaticImages{M, T} <: AbstractVector{T}
    v::NTuple{M, T}
end
Base.size(::StaticImages{M}) where {M} = (M,)
Base.@propagate_inbounds Base.getindex(s::StaticImages, i::Int) = s.v[i]

# (Callable structs are outside DispatchDoctor's reach; JET's @test_opt checks this one.)
@inline (m::MulBy{Nimber{K, T}, T})(x::Nimber{K, T}) where {K, T <: Unsigned} =
    unsafe_frombits(Nimber{K, T}, apply_bytes(m.tab, bits(x)))

# --- Constant times array in broadcasts (J8c, J8d) --------------------------------------------

# Operands of another type (elements of a lazy argument in a broadcast) get the scalar product.
@inline (m::MulBy)(x) = m.c * x

"""
    ConstMulStrategy

How a broadcast of a constant times an array (`c .* a`, `a .* c`, and `c * a`, which Base
computes by broadcasting) is computed, chosen per field by [`default_constmul_strategy`](@ref):

- [`ScalarConstMul`](@ref): the scalar `*`, element by element;
- [`LazyMulBy`](@ref): the byte tables of [`MulBy`](@ref) inside the fused broadcast, for long
  arrays.

Measurements and the choice of the defaults: `src/bulk.md`, "Constant times array".
"""
abstract type ConstMulStrategy end

"""
    ScalarConstMul()

`c .* a` multiplies element by element with the scalar `*`; fused expressions stay fused. Base's
behaviour; no field uses it by default any more (it was the default up to ``\\mathbb{F}_{256}``
until J8e, `src/bulk.md`).
"""
struct ScalarConstMul <: ConstMulStrategy end

"""
    LazyMulBy(threshold::Int)
    LazyMulBy(F) = LazyMulBy(mulby_threshold(F))

`c .* a` stays lazy and is decided when the (fused) broadcast is materialised: if the result has
at least `threshold` elements, every product with `c` is a lookup in the byte tables of
[`MulBy`](@ref)`(c)`, built once per broadcast; below, it is Base's broadcast of `*`
([`ScalarConstMul`](@ref)), so short arrays cost nothing extra. Fused expressions stay fused.
Implemented as a broadcast style (`ConstMulStyle`) that wraps the style of the arguments. The
default for every field, with the measured break-even [`mulby_threshold`](@ref).
"""
struct LazyMulBy <: ConstMulStrategy
    threshold::Int
end
LazyMulBy(::Type{N}) where {N <: Nimber} = LazyMulBy(mulby_threshold(N))

"""
    default_constmul_strategy(F)

The [`ConstMulStrategy`](@ref) behind `c .* a` and `a .* c` for a constant `c::F` and an array
(or lazy broadcast) `a`: [`LazyMulBy`](@ref)`(F)` for every field (`src/bulk.md`). To change it,
redefine it for a field, e.g.
`FiniteFieldArithmetic.default_constmul_strategy(::Type{GF2_64}) = ScalarConstMul()`.
"""
default_constmul_strategy(::Type{N}) where {N <: Nimber} = LazyMulBy(N)

"""
    mulby_threshold(F)

The array length from which [`LazyMulBy`](@ref) builds tables: the smallest measured length from
which it beat the scalar `*` in `c .* a`, `out .= c .* a` and `out .= c .* a .+ b`
(`bench/constmul.jl`; `src/bulk.md`, "Constant times array"): 96 at GF(2) … GF(256), 768 at
GF(2^16), 192 at GF(2^32), 256 at GF(2^64), 512 at GF(2^128).
"""
mulby_threshold(::Type{<:Nimber}) = 96
mulby_threshold(::Type{GF2_16}) = 768
mulby_threshold(::Type{GF2_32}) = 192
mulby_threshold(::Type{GF2_64}) = 256
mulby_threshold(::Type{GF2_128}) = 512

const ArrayOrLazy = Union{AbstractArray, Broadcast.Broadcasted}

"""
    const_mul_broadcasted(strategy::ConstMulStrategy, c, a)

What `Broadcast.broadcasted(*, c, a)` returns for a constant `c::Nimber` and an array or lazy
broadcast `a` under `strategy`: a lazy `Broadcasted`.
"""
const_mul_broadcasted(s::ConstMulStrategy, c::Nimber, a::ArrayOrLazy) =
    const_mul_broadcasted(s, c, a, (c, a))
const_mul_broadcasted(s::ConstMulStrategy, a::ArrayOrLazy, c::Nimber, args::Tuple) =
    const_mul_broadcasted(s, c, a, args)

# `args` is the argument tuple as written, (c, a) or (a, c); Base's broadcast keeps its order.
const_mul_broadcasted(::ScalarConstMul, ::Nimber, ::ArrayOrLazy, args::Tuple) =
    Broadcast.broadcasted(Broadcast.combine_styles(args...), *, args...)
const_mul_broadcasted(s::LazyMulBy, c::Nimber, a::ArrayOrLazy, args::Tuple) =
    Broadcast.Broadcasted{ConstMulStyle{typeof(Broadcast.combine_styles(args...))}}(
        ConstTimes(c, s.threshold), (a,))

# --- LazyMulBy as a broadcast style ---
#
# The decision between scalar products and tables must be one branch with two concrete call
# sites around the whole fused loop: a Union of two Broadcasted types returned by the hook (split
# by the compiler at materialize!) measured 1.74 ns per element at GF(2^16) fused, an explicit
# branch 1.50, Base's scalar loop 1.50 (src/bulk.md). The hook therefore returns a marker node,
# `ConstTimes`, in the style ConstMulStyle{S} (S: the style the arguments would have had). The
# style propagates to every enclosing node of the fused expression; `copyto!` and `copy` then
# compare the length of the result with the thresholds once and rewrite the tree to style S, with
# each marker replaced by Base's `*` broadcast or by a broadcast of MulBy(c).

# c ⊗ x, and the threshold of LazyMulBy. Callable, so an unrewritten tree is still correct.
struct ConstTimes{N <: Nimber} <: Function
    c::N
    threshold::Int
end
@inline (t::ConstTimes)(x) = t.c * x

struct ConstMulStyle{S <: Broadcast.BroadcastStyle} <: Broadcast.BroadcastStyle end

constmul_style(::Type{S}) where {S} = ConstMulStyle{S}()
Broadcast.BroadcastStyle(::ConstMulStyle{S}, o::Broadcast.BroadcastStyle) where {S} =
    constmul_style(typeof(Broadcast.result_style(S(), o)))
Broadcast.BroadcastStyle(o::Broadcast.BroadcastStyle, ::ConstMulStyle{S}) where {S} =
    constmul_style(typeof(Broadcast.result_style(o, S())))
Broadcast.BroadcastStyle(::ConstMulStyle{S}, ::ConstMulStyle{R}) where {S, R} =
    constmul_style(typeof(Broadcast.result_style(S(), R())))
# (Base's Unknown style would otherwise win against ours; keep ours, wrapping it.)
Broadcast.BroadcastStyle(::ConstMulStyle{S}, ::Broadcast.Unknown) where {S} = ConstMulStyle{S}()
Broadcast.BroadcastStyle(::Broadcast.Unknown, ::ConstMulStyle{S}) where {S} = ConstMulStyle{S}()

# The rewrite: Val(true) uses tables at every marker, Val(false) Base's `*`.
unwrap_constmul(x, ::Val) = x
unwrap_constmul(bc::Broadcast.Broadcasted{ConstMulStyle{S}}, v::Val) where {S} =
    Broadcast.Broadcasted{S}(unwrap_node(bc.f, unwrap_args(bc.args, v), v)..., bc.axes)
unwrap_args(args::Tuple, v::Val) = map(x -> unwrap_constmul(x, v), args)
# (function, arguments) of a rewritten node
unwrap_node(f, args::Tuple, ::Val) = (f, args)
unwrap_node(t::ConstTimes, args::Tuple, ::Val{true}) = (MulBy(t.c), args)
unwrap_node(t::ConstTimes, args::Tuple, ::Val{false}) = (*, (t.c, args...))

# The largest threshold of the markers in a tree (tables are used only if the result reaches it).
max_threshold(_) = 0
max_threshold(bc::Broadcast.Broadcasted) = max(node_threshold(bc.f), max_threshold_args(bc.args))
max_threshold_args(args::Tuple) = mapreduce(max_threshold, max, args; init = 0)
node_threshold(_) = 0
node_threshold(t::ConstTimes) = t.threshold

# Each branch behind a function barrier: with both loops inlined into one method, the fused loop
# at GF(2^16) ran at 1.63–1.65 ns per element against 1.47–1.49 for Base's (bench, src/bulk.md);
# with the barrier, equal from 1024 elements, about 9 ns more per call below.
@noinline unwrapped_copyto!(dest::AbstractArray, bc::Broadcast.Broadcasted, v::Val) =
    copyto!(dest, Broadcast.instantiate(unwrap_constmul(bc, v)))
@noinline unwrapped_copy(bc::Broadcast.Broadcasted, v::Val) =
    copy(Broadcast.instantiate(unwrap_constmul(bc, v)))

Base.copyto!(dest::AbstractArray, bc::Broadcast.Broadcasted{<:ConstMulStyle}) =
    length(dest) >= max_threshold(bc) ? unwrapped_copyto!(dest, bc, Val(true)) :
                                        unwrapped_copyto!(dest, bc, Val(false))
Base.copy(bc::Broadcast.Broadcasted{<:ConstMulStyle}) =
    length(bc) >= max_threshold(bc) ? unwrapped_copy(bc, Val(true)) : unwrapped_copy(bc, Val(false))

# The hooks. Only arrays whose element type is that of the constant, and lazy arguments (whose
# element type is not known here; MulBy falls back to `*` for other operands).
Broadcast.broadcasted(::typeof(*), c::N, a::AbstractArray{N}) where {N <: Nimber} =
    const_mul_broadcasted(default_constmul_strategy(N), c, a, (c, a))
Broadcast.broadcasted(::typeof(*), a::AbstractArray{N}, c::N) where {N <: Nimber} =
    const_mul_broadcasted(default_constmul_strategy(N), c, a, (a, c))
Broadcast.broadcasted(::typeof(*), c::N, a::Broadcast.Broadcasted) where {N <: Nimber} =
    const_mul_broadcasted(default_constmul_strategy(N), c, a, (c, a))
Broadcast.broadcasted(::typeof(*), a::Broadcast.Broadcasted, c::N) where {N <: Nimber} =
    const_mul_broadcasted(default_constmul_strategy(N), c, a, (a, c))

# --- Matrix products --------------------------------------------------------------------------

"""
    MatMulStrategy

How a matrix product over a nimber field is computed by [`matmul!`](@ref) (and hence by `*` and
`LinearAlgebra.mul!`, through [`default_matmul_strategy`](@ref)). All strategies give the same
result:

- [`ScalarMatMul`](@ref): the textbook loops with a scalar [`MulStrategy`](@ref).
- [`TableMatMul`](@ref): multiplication-by-constant tables ([`MulBy`](@ref)) per entry of `B`.
- [`CLMulMatMul`](@ref): carry-less products in the polynomial basis with one reduction per
  entry of the result.
"""
abstract type MatMulStrategy end

"""
    ScalarMatMul(strategy::MulStrategy)

The product ``C_{ij} = \\sum_k A_{ik} B_{kj}`` by the textbook loops (column by column, `B[k, j]`
times the column `A[:, k]`), each product with `strategy`.
"""
struct ScalarMatMul{S <: MulStrategy} <: MatMulStrategy
    strategy::S
end

"""
    TableMatMul()

The product column by column: ``C_{:,j} = \\sum_k A_{:,k} \\otimes B_{kj}``, each term through the
byte tables of multiplication by ``B_{kj}`` ([`MulBy`](@ref)), built once per entry of `B` into
a reused buffer and applied to the whole column of `A`. Up to ``\\mathbb{F}_{256}`` the table is
a column of the GF(256) product table and costs nothing to set up; above, the set-up (``n``
products) is paid per entry of `B` and only amortised over tall `A` (`src/bulk.md`).
"""
struct TableMatMul <: MatMulStrategy end

"""
    CLMulMatMul()

The product in the polynomial basis of [`CLMul`](@ref) (``\\mathbb{F}_{2^n}``, ``n \\ge 16``): both
matrices are converted once, every entry of the result is a sum of carry-less products, which is
reduced modulo ``f`` and converted back once (reduction is linear, so it commutes with the sum).
Needs `PCLMULQDQ` ([`has_pclmul`](@ref)); without it it computes as
`ScalarMatMul(MulTable())`, and up to ``\\mathbb{F}_{256}`` as [`TableMatMul`](@ref).
"""
struct CLMulMatMul <: MatMulStrategy end

"""
    default_matmul_strategy(F)

The [`MatMulStrategy`](@ref) behind `*` and `LinearAlgebra.mul!` for matrices over `F`, chosen by
the benchmarks in `bench/` (`src/bulk.md`, "Measurements").
"""
default_matmul_strategy(::Type{<:Nimber}) = TableMatMul()
default_matmul_strategy(::Type{<:Union{GF2_16, GF2_32, GF2_64, GF2_128}}) = CLMulMatMul()

"""
    matmul!(strategy::MatMulStrategy, C, A, B) -> C

Store the matrix product `A * B` in `C` (overwriting it), computed with `strategy`. `B` and `C`
may be vectors. `C` must not share memory with `A` or `B`.
"""
function matmul!(s::MatMulStrategy, C::StridedVecOrMat{N}, A::StridedMatrix{N},
                 B::StridedVecOrMat{N}) where {N <: Nimber}
    check_matmul(C, A, B)
    fill!(C, zero(N))
    return matmul_add!(s, C, A, B)
end

function check_matmul(C, A, B)
    m, p = size(A)
    (size(B, 1) == p && size(C, 1) == m && size(C, 2) == size(B, 2) && ndims(C) == ndims(B)) ||
        throw(DimensionMismatch("matmul!: C $(size(C)) ≠ A $(size(A)) × B $(size(B))"))
    (Base.mightalias(C, A) || Base.mightalias(C, B)) &&
        throw(ArgumentError("matmul!: output C must not share memory with A or B"))
    return nothing
end

@stable default_mode = "error" begin
    # C += A B, ScalarMatMul.
    function matmul_add!(sm::ScalarMatMul, C::StridedVecOrMat{N}, A::StridedMatrix{N},
                         B::StridedVecOrMat{N}) where {N <: Nimber}
        s = sm.strategy
        @inbounds for j in axes(B, 2), k in axes(A, 2)
            b = B[k, j]
            iszero(b) && continue
            for i in axes(A, 1)
                C[i, j] += mul(s, A[i, k], b)
            end
        end
        return C
    end

    # C += A B, TableMatMul.
    function matmul_add!(::TableMatMul, C::StridedVecOrMat{N}, A::StridedMatrix{N},
                         B::StridedVecOrMat{N}) where {K, T, N <: Nimber{K, T}}
        tab = Matrix{T}(undef, 256, sizeof(T))
        @inbounds for j in axes(B, 2), k in axes(A, 2)
            b = B[k, j]
            iszero(b) && continue
            mulby_tables!(tab, b)
            for i in axes(A, 1)
                C[i, j] += unsafe_frombits(N, apply_bytes(tab, bits(A[i, k])))
            end
        end
        return C
    end

    function matmul_add!(::CLMulMatMul, C::StridedVecOrMat{N}, A::StridedMatrix{N},
                         B::StridedVecOrMat{N}) where {K, T, N <: Nimber{K, T}}
        K <= 3 && return matmul_add!(TableMatMul(), C, A, B)
        has_pclmul() || return matmul_add!(ScalarMatMul(MulTable()), C, A, B)
        m, p = size(A)
        q = size(B, 2)
        At = Matrix{T}(undef, p, m)           # A transposed, so both factors run along k
        Bp = Matrix{T}(undef, p, q)
        @inbounds for i in 1:m, k in 1:p
            At[k, i] = to_poly(bits(A[i, k]))
        end
        @inbounds for j in 1:q, k in 1:p
            Bp[k, j] = to_poly(bits(B[k, j]))
        end
        @inbounds for j in 1:q, i in 1:m
            r = poly_dot(view(At, :, i), view(Bp, :, j))
            C[i, j] += unsafe_frombits(N, from_poly(r))
        end
        return C
    end
end

# The carry-less product of the hardware when the caller has checked `has_pclmul()`; the branch
# of `clmul64` stays out of the inner loops.
@static if Sys.ARCH === :x86_64
    @inline clmul64_unchecked(a::UInt64, b::UInt64) = clmul64_hw(a, b)
else
    @inline clmul64_unchecked(a::UInt64, b::UInt64) = clmul64_sw(a, b)
end

@stable default_mode = "error" begin
    # Σ_k a_k b_k in F_2[x]/(f): the unreduced products are summed (XOR), then reduced once.
    @inline function poly_dot(a::AbstractVector{T}, b::AbstractVector{T}) where {T <: Union{UInt16, UInt32}}
        acc = zero(UInt64)                     # products have degree < 2n - 1 ≤ 63
        @inbounds for k in eachindex(a, b)
            acc ⊻= clmul64_unchecked(a[k] % UInt64, b[k] % UInt64) % UInt64
        end
        n = 8 * sizeof(T)
        return poly_reduce((acc >> n) % T, acc % T)
    end
    @inline function poly_dot(a::AbstractVector{UInt64}, b::AbstractVector{UInt64})
        acc = zero(UInt128)
        @inbounds for k in eachindex(a, b)
            acc ⊻= clmul64_unchecked(a[k], b[k])
        end
        return poly_reduce((acc >> 64) % UInt64, acc % UInt64)
    end
    # 128 bits: Karatsuba on 64-bit halves is linear in each product, so its three partial
    # products are summed separately and combined once.
    @inline function poly_dot(a::AbstractVector{UInt128}, b::AbstractVector{UInt128})
        l = h = m = zero(UInt128)
        @inbounds for k in eachindex(a, b)
            x, y = a[k], b[k]
            x1, x0 = (x >> 64) % UInt64, x % UInt64
            y1, y0 = (y >> 64) % UInt64, y % UInt64
            l ⊻= clmul64_unchecked(x0, y0)
            h ⊻= clmul64_unchecked(x1, y1)
            m ⊻= clmul64_unchecked(x0 ⊻ x1, y0 ⊻ y1)
        end
        m ⊻= l ⊻ h
        return poly_reduce(h ⊻ (m >> 64), l ⊻ (m << 64))
    end
end

# `*` and `mul!` for strided matrices and vectors over a nimber field: C = A B α + C β, with α, β
# Booleans or field elements (integers other than Booleans are not field elements here: an
# integer is a bit pattern, not the image of ℤ, see binaryfield.jl). Matrix–matrix and
# matrix–vector are separate methods, so that each is more specific than LinearAlgebra's.
const Scalar{N} = Union{Bool, N}

function LinearAlgebra.mul!(C::StridedMatrix{N}, A::StridedMatrix{N}, B::StridedMatrix{N},
                            α::Scalar{N}, β::Scalar{N}) where {N <: Nimber}
    return mul_add_scaled!(C, A, B, α, β)
end
function LinearAlgebra.mul!(C::StridedVector{N}, A::StridedMatrix{N}, B::StridedVector{N},
                            α::Scalar{N}, β::Scalar{N}) where {N <: Nimber}
    return mul_add_scaled!(C, A, B, α, β)
end

function mul_add_scaled!(C::StridedVecOrMat{N}, A, B, α, β) where {N <: Nimber}
    check_matmul(C, A, B)
    s = default_matmul_strategy(N)
    if iszero(β)
        fill!(C, zero(N))
    elseif !isone(β)
        C .*= β
    end
    iszero(α) && return C
    isone(α) && return matmul_add!(s, C, A, B)
    P = matmul!(s, similar(C), A, B)
    C .+= P .* α
    return C
end

# --- Bitsliced vectors ------------------------------------------------------------------------

# Four 64-bit words processed as one SIMD vector (256 bits; LLVM splits or merges it to the
# target's vector width). Only AND and XOR are needed.
const LANES = 4
const Vec4 = NTuple{LANES, VecElement{UInt64}}
struct Word4
    v::Vec4
end
for (f, op) in ((:⊻, "xor"), (:&, "and"))
    ir = "%r = $op <$LANES x i64> %0, %1\nret <$LANES x i64> %r"
    @eval @inline Base.$f(a::Word4, b::Word4) =
        Word4(Base.llvmcall($ir, Vec4, Tuple{Vec4, Vec4}, a.v, b.v))
end

"""
    BitslicedVector{N}(undef, n)
    BitslicedVector(a::AbstractVector{N})

A vector of ``n`` elements of the nimber field `N` (up to ``\\mathbb{F}_{2^{32}}``) stored
*bitsliced*: plane ``j`` holds bit ``j`` of every element, 64 elements per word. Element-wise
products are then a fixed circuit of AND and XOR on whole words (the recursion of
[`Recursive`](@ref), with ``\\mathbb{F}_2`` products as AND), run on four words at a time with
SIMD, so one circuit evaluation multiplies 256 pairs. Element-wise `.*` and `.+` of two
`BitslicedVector`s and [`mul_elementwise!`](@ref) use it; everything else (indexing, other
broadcasts) works element by element and is slow. Convert once with `BitslicedVector(a)` and back
with `Vector(v)`. Measurements in `src/bulk.md`.
"""
struct BitslicedVector{N <: Nimber} <: AbstractVector{N}
    planes::Matrix{UInt64}     # planes[w, j]: bit j - 1 of elements 64(w - 1) + 1 … 64w
    len::Int
    function BitslicedVector{N}(::UndefInitializer, n::Integer) where {N <: Nimber}
        degree(N) <= 32 || throw(ArgumentError(
            "BitslicedVector: GF(2^$(degree(N))) is not supported (at most GF(2^32))"))
        nw = LANES * cld(cld(n, 64), LANES)    # whole SIMD words; the padding stays zero
        return new{N}(zeros(UInt64, nw, degree(N)), n)
    end
end

Base.size(v::BitslicedVector) = (v.len,)
Base.IndexStyle(::Type{<:BitslicedVector}) = IndexLinear()
Base.similar(v::BitslicedVector{N}) where {N} = BitslicedVector{N}(undef, length(v))

Base.@propagate_inbounds function Base.getindex(v::BitslicedVector{N}, i::Int) where {N}
    @boundscheck checkbounds(v, i)
    T = storagetype(N)
    w, t = ((i - 1) >> 6) + 1, (i - 1) & 63
    x = zero(T)
    @inbounds for j in 1:degree(N)
        x |= (((v.planes[w, j] >> t) & one(UInt64)) % T) << (j - 1)
    end
    return unsafe_frombits(N, x)
end

Base.@propagate_inbounds function Base.setindex!(v::BitslicedVector{N}, x, i::Int) where {N}
    @boundscheck checkbounds(v, i)
    b = bits(convert(N, x))
    w, t = ((i - 1) >> 6) + 1, (i - 1) & 63
    @inbounds for j in 1:degree(N)
        bit = UInt64((b >> (j - 1)) & one(b))
        v.planes[w, j] = (v.planes[w, j] & ~(one(UInt64) << t)) | (bit << t)
    end
    return v
end

function BitslicedVector(a::AbstractVector{N}) where {N <: Nimber}
    v = BitslicedVector{N}(undef, length(a))
    pack!(v.planes, a)
    return v
end

function Base.Vector{N}(v::BitslicedVector{N}) where {N}
    a = Vector{N}(undef, length(v))
    unpack!(a, v.planes)
    return a
end
Base.Vector(v::BitslicedVector{N}) where {N} = Vector{N}(v)
Base.collect(v::BitslicedVector{N}) where {N} = Vector{N}(v)

function Base.copyto!(dest::BitslicedVector{N}, src::BitslicedVector{N}) where {N}
    length(dest) == length(src) || throw(DimensionMismatch("copyto!: lengths differ"))
    copyto!(dest.planes, src.planes)
    return dest
end

# Element-wise operations on two bitsliced vectors stay bitsliced (`materialize` returns a
# non-`Broadcasted` value unchanged, so `a .* b` is this result).
function Base.Broadcast.broadcasted(::typeof(*), a::BitslicedVector{N},
                                    b::BitslicedVector{N}) where {N}
    return mul_elementwise!(similar(a), a, b)
end
function Base.Broadcast.broadcasted(::typeof(+), a::BitslicedVector{N},
                                    b::BitslicedVector{N}) where {N}
    length(a) == length(b) || throw(DimensionMismatch("lengths $(length(a)) and $(length(b))"))
    c = similar(a)
    c.planes .= a.planes .⊻ b.planes
    return c
end

"""
    mul_elementwise!(c, a, b) -> c

The element-wise product `c .= a .* b`. For [`BitslicedVector`](@ref)s it runs the bitsliced
circuit on the planes; for other arrays it is `broadcast!(*, c, a, b)`.
"""
mul_elementwise!(c::AbstractArray, a::AbstractArray, b::AbstractArray) = broadcast!(*, c, a, b)

@stable default_mode = "error" begin
    function mul_elementwise!(c::BitslicedVector{N}, a::BitslicedVector{N},
                              b::BitslicedVector{N}) where {N}
        length(a) == length(b) == length(c) ||
            throw(DimensionMismatch("lengths $(length(c)), $(length(a)), $(length(b))"))
        bitsliced_mul!(c.planes, a.planes, b.planes, Val(degree(N)))
        return c
    end
end

@inline load4(P::Matrix{UInt64}, w::Int, j::Int) =
    @inbounds Word4(ntuple(l -> VecElement(P[w + l - 1, j]), Val(LANES)))
@inline function store4!(P::Matrix{UInt64}, x::Word4, w::Int, j::Int)
    for l in 1:LANES
        @inbounds P[w + l - 1, j] = x.v[l].value
    end
    return nothing
end

function bitsliced_mul!(C::Matrix{UInt64}, A::Matrix{UInt64}, B::Matrix{UInt64},
                        ::Val{n}) where {n}
    for w in 1:LANES:size(A, 1)
        a = ntuple(j -> load4(A, w, j), Val(n))
        b = ntuple(j -> load4(B, w, j), Val(n))
        r = bs_mul(a, b)
        for j in 1:n
            store4!(C, r[j], w, j)
        end
    end
    return C
end

# The recursion of `Recursive` on planes: an element of level K is a tuple of 2^K planes (lowest
# bit first), its halves are the two halves of the tuple, ⊕ is XOR of planes, and the level-0
# product is AND. The tuples are fully unrolled, so the whole circuit is straight-line code.
# One method per tuple length (levels 1–5), so that every signature is concrete.
const Planes{M} = NTuple{M, Word4}

@inline bs_half(c::Planes{1}) = c
@inline bs_mul(a::Planes{1}, b::Planes{1}) = (a[1] & b[1],)
for M in (2, 4, 8, 16, 32)
    h = M ÷ 2
    @eval begin
        @inline lo_half(t::Planes{$M}) = ntuple(i -> t[i], Val($h))
        @inline hi_half(t::Planes{$M}) = ntuple(i -> t[$h + i], Val($h))
        @inline plane_xor(a::Planes{$h}, b::Planes{$h}) = map(⊻, a, b)
        @inline function bs_half(c::Planes{$M})                 # c ⊗ 2^(2^k - 1), see nimber.jl
            c0, c1 = lo_half(c), hi_half(c)
            return (bs_half(bs_half(c1))..., bs_half(plane_xor(c0, c1))...)
        end
        @inline function bs_mul(a::Planes{$M}, b::Planes{$M})
            a0, a1, b0, b1 = lo_half(a), hi_half(a), lo_half(b), hi_half(b)
            c = bs_mul(a1, b1)
            d = bs_mul(a0, b0)
            e = bs_mul(plane_xor(a0, a1), plane_xor(b0, b1))
            return (plane_xor(d, bs_half(c))..., plane_xor(e, d)...)
        end
    end
end

# Packing: plane j of word w collects bit j of 64 consecutive elements. The inner loop runs over
# the 64 elements (variable shifts, vectorisable); elements past the end count as zero.
function pack!(P::Matrix{UInt64}, a::AbstractVector{N}) where {N <: Nimber}
    n = length(a)
    T = storagetype(N)
    buf = zeros(T, 64)
    for w in 1:cld(n, 64)
        base = 64 * (w - 1)
        m = min(64, n - base)
        @inbounds for t in 1:m
            buf[t] = bits(a[firstindex(a) + base + t - 1])
        end
        @inbounds for t in (m + 1):64
            buf[t] = zero(T)
        end
        @inbounds for j in 1:degree(N)
            acc = zero(UInt64)
            for t in 1:64
                acc |= (((buf[t] >> (j - 1)) & one(T)) % UInt64) << (t - 1)
            end
            P[w, j] = acc
        end
    end
    return P
end

function unpack!(a::AbstractVector{N}, P::Matrix{UInt64}) where {N <: Nimber}
    n = length(a)
    T = storagetype(N)
    buf = zeros(T, 64)
    for w in 1:cld(n, 64)
        fill!(buf, zero(T))
        @inbounds for j in 1:degree(N)
            p = P[w, j]
            for t in 1:64
                buf[t] |= (((p >> (t - 1)) & one(UInt64)) % T) << (j - 1)
            end
        end
        base = 64 * (w - 1)
        @inbounds for t in 1:min(64, n - base)
            a[firstindex(a) + base + t - 1] = unsafe_frombits(N, buf[t])
        end
    end
    return a
end

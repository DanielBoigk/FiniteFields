# The bitsliced matrix product (backlog J8b): a matrix over GF(2^n) as n bit matrices (one per
# bit of the entries), the product as the nimber recursion of `Recursive` evaluated on whole bit
# matrices (3^K products over GF(2) at level K), each GF(2) product by the method of the Four
# Russians. Theory and measurements: bitsliced_matmul.md beside this file.

"""
    BitslicedMatMul()

The matrix product on bitsliced matrices: both factors are split into ``n`` bit matrices (plane
``j`` holds bit ``j`` of every entry), and the product is the recursion of [`Recursive`](@ref)
evaluated on whole matrices, since every step of it is bilinear or linear: ``3^K`` matrix
products over ``\\mathbb{F}_2`` at ``\\mathbb{F}_{2^{2^K}}``, each by the method of the Four
Russians (byte tables of sums of rows). Converts once per product. Fast for large matrices over
small fields; measurements in `src/bitsliced_matmul.md`.
"""
struct BitslicedMatMul <: MatMulStrategy end

# Layouts. A *left* bit-matrix array `L::Array{UInt64, 3}` of size (nw, rows, planes): `L[w, i, j]`
# holds bit j - 1 of the entries in row i, columns 64(w - 1) + 1 … 64w (column c at bit c - 1 of
# its word). A *right* array `R::Array{Word4, 3}` of size (nc, rows, planes): the same with 256
# columns per `Word4`. A = left (its columns are the inner dimension), B and C = right. Padding
# bits are zero; B has 64 nw(A) rows, the rows beyond the inner dimension zero.

const ZERO4 = Word4(ntuple(_ -> VecElement(zero(UInt64)), Val(LANES)))

@stable default_mode = "error" begin
    function matmul_add!(::BitslicedMatMul, C::StridedVecOrMat{N}, A::StridedMatrix{N},
                         B::StridedVecOrMat{N}) where {N <: Nimber}
        m, p = size(A)
        q = size(B, 2)
        (m == 0 || p == 0 || q == 0) && return C
        n = Val(degree(N))
        nw = cld(p, 64)
        L = left_planes(A, n)
        R = right_planes(B, 64 * nw, n)
        P = Array{Word4, 3}(undef, size(R, 1), m, degree(N))
        T = Vector{Word4}(undef, 8 * 256)
        bitsliced_matmul!(P, 0, L, 0, R, 0, degree(N), T)
        add_right_planes!(C, P, n)
        return C
    end

    # The planes of a matrix, rows packed: left layout.
    function left_planes(A::AbstractVecOrMat{N}, ::Val{n}) where {N <: Nimber, n}
        m, p = size(A, 1), size(A, 2)
        nw = cld(p, 64)
        L = Array{UInt64, 3}(undef, nw, m, n)
        buf = zeros(N, 64)
        for w in 1:nw, i in 1:m
            row_block!(buf, A, i, 64 * (w - 1))
            block_planes!(view(L, w, i, :), buf, 1, Val(n))
        end
        return L
    end

    # The planes of a matrix in the right layout, with `rows` ≥ size(B, 1) rows (zeros below).
    function right_planes(B::AbstractVecOrMat{N}, rows::Int, ::Val{n}) where {N <: Nimber, n}
        p, q = size(B, 1), size(B, 2)
        nc = cld(q, 64 * LANES)
        U = zeros(UInt64, LANES * nc, rows, n)
        buf = zeros(N, 64)
        for w in 1:cld(q, 64), k in 1:p
            row_block!(buf, B, k, 64 * (w - 1))
            block_planes!(view(U, w, k, :), buf, 1, Val(n))
        end
        R = Array{Word4, 3}(undef, nc, rows, n)
        for j in 1:n, k in 1:rows, c in 1:nc
            @inbounds R[c, k, j] = Word4(ntuple(l -> VecElement(U[LANES * (c - 1) + l, k, j]),
                                                Val(LANES)))
        end
        return R
    end

    # C += the matrix whose planes are P (right layout).
    function add_right_planes!(C::AbstractVecOrMat{N}, P::Array{Word4, 3}, ::Val{n}) where {N, n}
        m, q = size(C, 1), size(C, 2)
        buf = Vector{N}(undef, 64)
        y = Vector{UInt64}(undef, n)
        for c in axes(P, 1), l in 1:LANES, i in 1:m
            c0 = 64 * (LANES * (c - 1) + l - 1)
            c0 >= q && break
            for j in 1:n
                @inbounds y[j] = P[c, i, j].v[l].value
            end
            block_from_planes!(buf, 1, y, Val(n))
            for t in 1:min(64, q - c0)
                @inbounds C[i, c0 + t] += buf[t]
            end
        end
        return C
    end

    # out[1:64] = X[i, c0 + 1 : c0 + 64], zeros beyond the last column.
    @inline function row_block!(buf::Vector{N}, X::AbstractVecOrMat{N}, i::Int, c0::Int) where {N}
        len = clamp(size(X, 2) - c0, 0, 64)
        for t in 1:len
            @inbounds buf[t] = X[i, c0 + t]
        end
        for t in (len + 1):64
            @inbounds buf[t] = zero(N)
        end
        return buf
    end

    # P[:, :, po + 1 : po + np] = the product of the level-log2(np) bit-matrix tuples
    # L[:, :, lo + 1 : lo + np] and R[:, :, ro + 1 : ro + np]: the recursion of `Recursive`
    # (`bs_mul` in bulk.jl) with matrices for planes. Temporaries are allocated per level.
    function bitsliced_matmul!(P::Array{Word4, 3}, po::Int, L::Array{UInt64, 3}, lo::Int,
                               R::Array{Word4, 3}, ro::Int, np::Int, T::Vector{Word4})
        if np == 1
            fill!(view(P, :, :, po + 1), ZERO4)
            gf2_matmul_add!(P, po + 1, L, lo + 1, R, ro + 1, T)
            return P
        end
        h = np >> 1
        nc, m = size(P, 1), size(P, 2)
        bitsliced_matmul!(P, po, L, lo, R, ro, h, T)                          # d = a0 b0
        Cc = Array{Word4, 3}(undef, nc, m, h)
        bitsliced_matmul!(Cc, 0, L, lo + h, R, ro + h, h, T)                  # c = a1 b1
        SL = Array{UInt64, 3}(undef, size(L, 1), size(L, 2), h)
        SR = Array{Word4, 3}(undef, size(R, 1), size(R, 2), h)
        @inbounds for j in 1:h
            for x in axes(L, 2), w in axes(L, 1)
                SL[w, x, j] = L[w, x, lo + j] ⊻ L[w, x, lo + h + j]
            end
            for x in axes(R, 2), w in axes(R, 1)
                SR[w, x, j] = R[w, x, ro + j] ⊻ R[w, x, ro + h + j]
            end
        end
        bitsliced_matmul!(P, po + h, SL, 0, SR, 0, h, T)                       # e
        masks = HALF_MASKS[trailing_zeros(h) + 1]
        @inbounds for j in 1:h
            for x in 1:m, w in 1:nc
                P[w, x, po + h + j] ⊻= P[w, x, po + j]                         # e ⊕ d
            end
            s = masks[j]                                                       # d ⊕ half(c)
            while s != 0
                k = trailing_zeros(s) + 1
                s &= s - one(s)
                for x in 1:m, w in 1:nc
                    P[w, x, po + j] ⊻= Cc[w, x, k]
                end
            end
        end
        return P
    end

    # P[:, :, jp] ⊻= L[:, :, jl] R[:, :, jr] over GF(2), by the method of the Four Russians: per
    # 64 rows of R (one word of L) eight tables of the 256 sums of 8 rows each, then per row of
    # L eight lookups. Tables are built by doubling (entry 2^b + s = entry s ⊕ row b), so the
    # XORs of one step are independent.
    function gf2_matmul_add!(P::Array{Word4, 3}, jp::Int, L::Array{UInt64, 3}, jl::Int,
                             R::Array{Word4, 3}, jr::Int, T::Vector{Word4})
        m = size(P, 2)
        @inbounds for c in axes(P, 1), wl in axes(L, 1)
            g = 64 * (wl - 1)
            for t in 0:7
                base = 256 * t
                T[base + 1] = ZERO4
                for b in 0:7
                    r = R[c, g + 8 * t + b + 1, jr]
                    for s in 0:((1 << b) - 1)
                        T[base + (1 << b) + s + 1] = T[base + s + 1] ⊻ r
                    end
                end
            end
            for i in 1:m
                a = L[wl, i, jl]
                acc = P[c, i, jp]
                for t in 0:7
                    acc = acc ⊻ T[256 * t + ((a >> (8 * t)) & 0xff) % Int + 1]
                end
                P[c, i, jp] = acc
            end
        end
        return P
    end
end

# The linear map half of `bs_half` (bulk.jl) at level L as plane masks: plane j of half(c) is the
# XOR of the planes k of c with bit k - 1 of mask j set. Same recursion as `bs_half`, evaluated
# on masks instead of planes: half(c) = (half(half(c1)), half(c0 ⊕ c1)), half = id at level 0.
function half_masks(L::Int)
    L == 0 && return [UInt64(1)]
    h = 1 << (L - 1)
    sub = half_masks(L - 1)
    apply(ms, v) = [foldl(⊻, (v[k] for k in 1:length(v) if (ms[j] >> (k - 1)) & 1 == 1);
                          init = zero(UInt64)) for j in eachindex(ms)]
    c0 = [UInt64(1) << (k - 1) for k in 1:h]
    c1 = [UInt64(1) << (h + k - 1) for k in 1:h]
    return vcat(apply(sub, apply(sub, c1)), apply(sub, c0 .⊻ c1))
end
const HALF_MASKS = [half_masks(L) for L in 0:6]      # levels up to 6 (planes up to 64 per half)

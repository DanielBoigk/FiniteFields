module FiniteFields

export Fp

function _isprime(n::Integer)
    n < 2 && return false
    iseven(n) && return n == 2
    d = 3
    while d <= n ÷ d
        n % d == 0 && return false
        d += 2
    end
    true
end

"""
    Fp{P}(n)

The residue class of `n` in the prime field with `P` elements.
Construction checks that `P` is prime. This starter implementation is intended for small fields.
"""
struct Fp{P} <: Number
    value::Int

    function Fp{P}(n::Integer) where {P}
        (P isa Integer && 2 <= P <= typemax(Int)) ||
            throw(ArgumentError("field characteristic must be an integer between 2 and typemax(Int)"))
        _isprime(P) || throw(ArgumentError("field characteristic must be prime"))
        new{P}(Int(mod(big(n), big(P))))
    end
end

Base.zero(::Type{Fp{P}}) where {P} = Fp{P}(0)
Base.one(::Type{Fp{P}}) where {P} = Fp{P}(1)
Base.iszero(x::Fp) = iszero(x.value)
Base.:(==)(x::Fp{P}, y::Fp{P}) where {P} = x.value == y.value
Base.:+(x::Fp{P}, y::Fp{P}) where {P} = Fp{P}(big(x.value) + y.value)
Base.:-(x::Fp{P}, y::Fp{P}) where {P} = Fp{P}(big(x.value) - y.value)
Base.:*(x::Fp{P}, y::Fp{P}) where {P} = Fp{P}(big(x.value) * y.value)
Base.:-(x::Fp{P}) where {P} = Fp{P}(-x.value)

function Base.inv(x::Fp{P}) where {P}
    iszero(x) && throw(DivideError())
    Fp{P}(invmod(x.value, P))
end

Base.:/(x::Fp{P}, y::Fp{P}) where {P} = x * inv(y)
Base.show(io::IO, x::Fp{P}) where {P} = print(io, x.value, " (mod ", P, ")")

end

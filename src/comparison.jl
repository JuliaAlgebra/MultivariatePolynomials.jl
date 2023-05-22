export isapproxzero

export LexOrder, InverseLexOrder, Reverse, Graded

Base.iszero(v::AbstractVariable) = false
Base.iszero(m::AbstractMonomial) = false
Base.iszero(t::AbstractTerm) = iszero(coefficient(t))
Base.iszero(t::AbstractPolynomial) = iszero(nterms(t))

Base.isone(v::AbstractVariable) = false
Base.isone(m::AbstractMonomial) = isconstant(m)
Base.isone(t::AbstractTerm) = isone(coefficient(t)) && isconstant(monomial(t))
function Base.isone(p::AbstractPolynomial)
    return isone(nterms(p)) && isone(first(terms(p)))
end

# See https://github.com/blegat/MultivariatePolynomials.jl/issues/22
# avoids the call to be transfered to left_constant_eq
Base.:(==)(α::Nothing, x::APL) = false
Base.:(==)(x::APL, α::Nothing) = false
Base.:(==)(α::Dict, x::APL) = false
Base.:(==)(x::APL, α::Dict) = false
Base.:(==)(α::Nothing, x::RationalPoly) = false
Base.:(==)(x::RationalPoly, α::Nothing) = false
Base.:(==)(α::Dict, x::RationalPoly) = false
Base.:(==)(x::RationalPoly, α::Dict) = false

function right_term_eq(p::AbstractPolynomial, t)
    if iszero(p)
        iszero(t)
    else
        # terms/nterms ignore zero terms
        nterms(p) == 1 && leading_term(p) == t
    end
end
right_term_eq(p::APL, t) = right_term_eq(polynomial(p), t)

left_constant_eq(α, v::AbstractVariable) = false
right_constant_eq(v::AbstractVariable, α) = false
function _term_constant_eq(t::AbstractTermLike, α)
    if iszero(t)
        iszero(α)
    else
        α == coefficient(t) && isconstant(t)
    end
end
left_constant_eq(α, t::AbstractTermLike) = _term_constant_eq(t, α)
right_constant_eq(t::AbstractTermLike, α) = _term_constant_eq(t, α)
left_constant_eq(α, p::APL) = right_term_eq(p, α)
right_constant_eq(p::APL, α) = right_term_eq(p, α)

function Base.:(==)(mono::AbstractMonomial, v::AbstractVariable)
    return isone(degree(mono)) && variable(mono) == v
end
function Base.:(==)(v::AbstractVariable, mono::AbstractMonomial)
    return isone(degree(mono)) && v == variable(mono)
end
function Base.:(==)(t::AbstractTerm, mono::AbstractMonomialLike)
    return isone(coefficient(t)) && monomial(t) == mono
end
function Base.:(==)(mono::AbstractMonomialLike, t::AbstractTerm)
    return isone(coefficient(t)) && mono == monomial(t)
end

function Base.:(==)(t1::AbstractTerm, t2::AbstractTerm)
    c1 = coefficient(t1)
    c2 = coefficient(t2)
    if iszero(c1)
        iszero(c2)
    else
        c1 == c2 && monomial(t1) == monomial(t2)
    end
end
Base.:(==)(p::AbstractPolynomial, t::AbstractTerm) = right_term_eq(p, t)
Base.:(==)(t::AbstractTerm, p::AbstractPolynomial) = right_term_eq(p, t)

function compare_terms(p1::AbstractPolynomial, p2::AbstractPolynomial, isz, op)
    i1 = 1
    i2 = 1
    t1 = terms(p1)
    t2 = terms(p2)
    while true
        while i1 <= length(t1) && isz(coefficient(t1[i1]))
            i1 += 1
        end
        while i2 <= length(t2) && isz(coefficient(t2[i2]))
            i2 += 1
        end
        if i1 > length(t1) && i2 > length(t2)
            return true
        end
        if i1 > length(t1) || i2 > length(t2)
            return false
        end
        if !op(t1[i1], t2[i2])
            return false
        end
        i1 += 1
        i2 += 1
    end
end

# Can there be zero term in TypedPolynomials ?
#function (==)(p1::AbstractPolynomial, p2::AbstractPolynomial)
#    nterms(p1) != nterms(p2) && return false
#    for (t1, t2) in zip(terms(p1), terms(p2))
#        @assert !iszero(t1) && !iszero(t2) # There should be no zero term
#        if t1 != t2
#            return false
#        end
#    end
#    return true
#end
function Base.:(==)(p1::AbstractPolynomial, p2::AbstractPolynomial)
    return compare_terms(p1, p2, iszero, ==)
end

Base.:(==)(p::RationalPoly, q::RationalPoly) = p.num * q.den == q.num * p.den
# Solve ambiguity with (::PolyType, ::Any)
Base.:(==)(p::APL, q::RationalPoly) = p * q.den == q.num
Base.:(==)(q::RationalPoly, p::APL) = p == q
Base.:(==)(α, q::RationalPoly) = α * q.den == q.num
Base.:(==)(q::RationalPoly, α) = α == q

# α could be a JuMP affine expression
isapproxzero(α; ztol::Real = 0.0) = false
function isapproxzero(α::Number; ztol::Real = Base.rtoldefault(α, α, 0))
    return abs(α) <= ztol
end

isapproxzero(m::AbstractMonomialLike; kwargs...) = false
function isapproxzero(t::AbstractTermLike; kwargs...)
    return isapproxzero(coefficient(t); kwargs...)
end
function isapproxzero(p::APL; kwargs...)
    return all(term -> isapproxzero(term; kwargs...), terms(p))
end
isapproxzero(p::RationalPoly; kwargs...) = isapproxzero(p.num; kwargs...)

function Base.isapprox(t1::AbstractTerm, t2::AbstractTerm; kwargs...)
    return isapprox(coefficient(t1), coefficient(t2); kwargs...) &&
           monomial(t1) == monomial(t2)
end
function Base.isapprox(
    p1::AbstractPolynomial{S},
    p2::AbstractPolynomial{T};
    atol = 0,
    ztol::Real = iszero(atol) ? Base.rtoldefault(S, T, 0) : atol,
    kwargs...,
) where {S,T}
    return compare_terms(
        p1,
        p2,
        t -> isapproxzero(t; ztol = ztol),
        (x, y) -> isapprox(x, y; atol = atol, kwargs...),
    )
end

function Base.isapprox(p::RationalPoly, q::RationalPoly; kwargs...)
    return isapprox(p.num * q.den, q.num * p.den; kwargs...)
end
function Base.isapprox(p::RationalPoly, q::APL; kwargs...)
    return isapprox(p.num, q * p.den; kwargs...)
end
function Base.isapprox(p::APL, q::RationalPoly; kwargs...)
    return isapprox(p * q.den, q.num; kwargs...)
end
function Base.isapprox(q::RationalPoly{C}, α; kwargs...) where {C}
    return isapprox(q, constant_term(α, q.den); kwargs...)
end
function Base.isapprox(α, q::RationalPoly{C}; kwargs...) where {C}
    return isapprox(constant_term(α, q.den), q; kwargs...)
end

"""
    abstract type AbstractMonomialOrdering end

Abstract type for monomial ordering as defined in [CLO13, Definition 2.2.1, p. 55]

[CLO13] Cox, D., Little, J., & OShea, D.
*Ideals, varieties, and algorithms: an introduction to computational algebraic geometry and commutative algebra*.
Springer Science & Business Media, **2013**.
"""
abstract type AbstractMonomialOrdering end

"""
    compare(a, b, order::Type{<:AbstractMonomialOrdering})

Returns a negative number if `a < b`, a positive number if `a > b` and zero if `a == b`.
The comparison is done according to `order`.
"""
function compare end

"""
    struct LexOrder <: AbstractMonomialOrdering end

Lexicographic (Lex for short) Order often abbreviated as *lex* order as defined in [CLO13, Definition 2.2.3, p. 56]

The [`Graded`](@ref) version is often abbreviated as *grlex* order and is defined in [CLO13, Definition 2.2.5, p. 58]

[CLO13] Cox, D., Little, J., & OShea, D.
*Ideals, varieties, and algorithms: an introduction to computational algebraic geometry and commutative algebra*.
Springer Science & Business Media, **2013**.
"""
struct LexOrder <: AbstractMonomialOrdering end

function compare(exp1::AbstractVector{T}, exp2::AbstractVector{T}, ::Type{LexOrder}) where {T}
    if eachindex(exp1) != eachindex(exp2)
        throw(ArgumentError("Cannot compare exponent vectors `$exp1` and `$exp2` of different indices."))
    end
    @inbounds for i in eachindex(exp1)
        Δ = exp1[i] - exp2[i]
        if !iszero(Δ)
            return Δ
        end
    end
    return zero(T)
end

"""
    struct InverseLexOrder <: AbstractMonomialOrdering end

Inverse Lex Order defined in [CLO13, Exercise 2.2.6, p. 61] where it is abbreviated as *invlex*.
It corresponds to [`LexOrder`](@ref) but with the variables in reverse order.

The [`Graded`](@ref) version can be abbreviated as *grinvlex* order.
It is defined in [BDD13, Definition 2.1] where it is called *Graded xel order*.

[CLO13] Cox, D., Little, J., & OShea, D.
*Ideals, varieties, and algorithms: an introduction to computational algebraic geometry and commutative algebra*.
Springer Science & Business Media, **2013**.
[BDD13] Batselier, K., Dreesen, P., & De Moor, B.
*The geometry of multivariate polynomial division and elimination*.
SIAM Journal on Matrix Analysis and Applications, 34(1), 102-125, *2013*.
"""
struct InverseLexOrder <: AbstractMonomialOrdering end

function compare(exp1::AbstractVector{T}, exp2::AbstractVector{T}, ::Type{InverseLexOrder}) where {T}
    if eachindex(exp1) != eachindex(exp2)
        throw(ArgumentError("Cannot compare exponent vectors `$exp1` and `$exp2` of different indices."))
    end
    @inbounds for i in Iterators.Reverse(eachindex(exp1))
        Δ = exp1[i] - exp2[i]
        if !iszero(Δ)
            return Δ
        end
    end
    return zero(T)
end

"""
    struct Graded{O<:AbstractMonomialOrdering} <: AbstractMonomialOrdering
        same_degree_ordering::O
    end

Monomial ordering defined by:
* `degree(a) == degree(b)` then the ordering is determined by `same_degree_ordering`,
* otherwise, it is the ordering between the integers `degree(a)` and `degree(b)`.
"""
struct Graded{O<:AbstractMonomialOrdering} <: AbstractMonomialOrdering
    same_degree_ordering::O
end

_deg(exponents) = sum(exponents)
_deg(mono::AbstractMonomial) = degree(mono)

function compare(a, b, ::Type{Graded{O}}) where {O}
    deg_a = _deg(a)
    deg_b = _deg(b)
    if deg_a == deg_b
        return compare(a, b, O)
    else
        return deg_a - deg_b
    end
end

"""
    struct Reverse{O<:AbstractMonomialOrdering} <: AbstractMonomialOrdering
        reverse_order::O
    end

Monomial ordering defined by
`compare(a, b, ::Type{Reverse{O}}) where {O} = compare(b, a, O)`.

Reverse Lex Order defined in [CLO13, Exercise 2.2.9, p. 61] where it is abbreviated as *rinvlex*.
can be obtained as `Reverse(InverseLexOrder())`.

The Graded Reverse Lex Order often abbreviated as *grevlex* order defined in [CLO13, Definition 2.2.6, p. 58]
can be obtained as `Graded(Reverse(InverseLexOrder()))`.

[CLO13] Cox, D., Little, J., & OShea, D.
*Ideals, varieties, and algorithms: an introduction to computational algebraic geometry and commutative algebra*.
Springer Science & Business Media, **2013**.
"""
struct Reverse{O<:AbstractMonomialOrdering} <: AbstractMonomialOrdering
    reverse_ordering::O
end

compare(a, b, ::Type{Reverse{O}}) where {O} = compare(b, a, O)

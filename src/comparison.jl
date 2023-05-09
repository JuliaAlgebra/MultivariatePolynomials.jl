export isapproxzero

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
Base.:(==)(p1::AbstractPolynomial, p2::AbstractPolynomial) = compare_terms(p1, p2, iszero, ==)

Base.:(==)(p::RationalPoly, q::RationalPoly) = p.num*q.den == q.num*p.den
# Solve ambiguity with (::PolyType, ::Any)
Base.:(==)(p::APL, q::RationalPoly) = p*q.den == q.num
Base.:(==)(q::RationalPoly, p::APL) = p == q
Base.:(==)(α, q::RationalPoly) = α*q.den == q.num
Base.:(==)(q::RationalPoly, α) = α == q

# α could be a JuMP affine expression
isapproxzero(α; ztol::Real=0.) = false
function isapproxzero(α::Number; ztol::Real=Base.rtoldefault(α, α, 0))
    abs(α) <= ztol
end

isapproxzero(m::AbstractMonomialLike; kwargs...) = false
isapproxzero(t::AbstractTermLike; kwargs...) = isapproxzero(coefficient(t); kwargs...)
isapproxzero(p::APL; kwargs...) = all(term -> isapproxzero(term; kwargs...), terms(p))
isapproxzero(p::RationalPoly; kwargs...) = isapproxzero(p.num; kwargs...)

Base.isapprox(t1::AbstractTerm, t2::AbstractTerm; kwargs...) = isapprox(coefficient(t1), coefficient(t2); kwargs...) && monomial(t1) == monomial(t2)
function Base.isapprox(p1::AbstractPolynomial{S}, p2::AbstractPolynomial{T};
                       atol=0, ztol::Real=iszero(atol) ? Base.rtoldefault(S, T, 0) : atol, kwargs...) where {S, T}
    return compare_terms(p1, p2, t -> isapproxzero(t; ztol=ztol),
                         (x, y) -> isapprox(x, y; atol=atol, kwargs...))
end

Base.isapprox(p::RationalPoly, q::RationalPoly; kwargs...) = isapprox(p.num*q.den, q.num*p.den; kwargs...)
Base.isapprox(p::RationalPoly, q::APL; kwargs...) = isapprox(p.num, q*p.den; kwargs...)
Base.isapprox(p::APL, q::RationalPoly; kwargs...) = isapprox(p*q.den, q.num; kwargs...)
Base.isapprox(q::RationalPoly{C}, α; kwargs...) where {C} = isapprox(q, constant_term(α, q.den); kwargs...)
Base.isapprox(α, q::RationalPoly{C}; kwargs...) where {C} = isapprox(constant_term(α, q.den), q; kwargs...)

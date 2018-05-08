export isapproxzero

Base.iszero(v::AbstractVariable) = false
Base.iszero(m::AbstractMonomial) = false
Base.iszero(t::AbstractTerm) = iszero(coefficient(t))
Base.iszero(t::AbstractPolynomial) = iszero(nterms(t))

# See https://github.com/blegat/MultivariatePolynomials.jl/issues/22
# avoids the call to be transfered to eqconstant
Base.:(==)(α::Nothing, x::APL) = false
Base.:(==)(x::APL, α::Nothing) = false
Base.:(==)(α::Dict, x::APL) = false
Base.:(==)(x::APL, α::Dict) = false
Base.:(==)(α::Nothing, x::RationalPoly) = false
Base.:(==)(x::RationalPoly, α::Nothing) = false
Base.:(==)(α::Dict, x::RationalPoly) = false
Base.:(==)(x::RationalPoly, α::Dict) = false

function polyeqterm(p::AbstractPolynomial, t)
    if iszero(p)
        iszero(t)
    else
        # terms/nterms ignore zero terms
        nterms(p) == 1 && leadingterm(p) == t
    end
end
polyeqterm(p::APL, t) = polyeqterm(polynomial(p), t)

eqconstant(α, v::AbstractVariable) = false
eqconstant(v::AbstractVariable, α) = false
function _termeqconstant(t::AbstractTermLike, α)
    if iszero(t)
        iszero(α)
    else
        α == coefficient(t) && isconstant(t)
    end
end
eqconstant(α, t::AbstractTermLike) = _termeqconstant(t, α)
eqconstant(t::AbstractTermLike, α) = _termeqconstant(t, α)
eqconstant(α, p::APL) = polyeqterm(p, α)
eqconstant(p::APL, α) = polyeqterm(p, α)

function Base.:(==)(t1::AbstractTerm, t2::AbstractTerm)
    c1 = coefficient(t1)
    c2 = coefficient(t2)
    if iszero(c1)
        iszero(c2)
    else
        c1 == c2 && monomial(t1) == monomial(t2)
    end
end
Base.:(==)(p::AbstractPolynomial, t::AbstractTerm) = polyeqterm(p, t)
Base.:(==)(t::AbstractTerm, p::AbstractPolynomial) = polyeqterm(p, t)

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
isapproxzero(p::APL; kwargs...) = all(isapproxzero.(terms(p); kwargs...))
isapproxzero(p::RationalPoly; kwargs...) = isapproxzero(p.num; kwargs...)

Base.isapprox(t1::AbstractTerm, t2::AbstractTerm; kwargs...) = isapprox(coefficient(t1), coefficient(t2); kwargs...) && monomial(t1) == monomial(t2)
Base.isapprox(p1::AbstractPolynomial, p2::AbstractPolynomial; ztol::Real=1e-6, kwargs...) = compare_terms(p1, p2, t -> isapproxzero(t; ztol=ztol), (x, y) -> isapprox(x, y; kwargs...))

Base.isapprox(p::RationalPoly, q::RationalPoly; kwargs...) = isapprox(p.num*q.den, q.num*p.den; kwargs...)
Base.isapprox(p::RationalPoly, q::APL; kwargs...) = isapprox(p.num, q*p.den; kwargs...)
Base.isapprox(p::APL, q::RationalPoly; kwargs...) = isapprox(p*q.den, q.num; kwargs...)
Base.isapprox(q::RationalPoly{C}, α; kwargs...) where {C} = isapprox(q, constantterm(α, q.den); kwargs...)
Base.isapprox(α, q::RationalPoly{C}; kwargs...) where {C} = isapprox(constantterm(α, q.den), q; kwargs...)

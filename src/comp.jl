Base.iszero(v::AbstractVariable) = false
Base.iszero(m::AbstractMonomial) = false
Base.iszero(t::AbstractTerm) = iszero(coefficient(t))
Base.iszero(p::MatPolynomial) = isempty(polynomial(p))

# See https://github.com/blegat/MultivariatePolynomials.jl/issues/22
# avoids the call to be transfered to eqconstant
(==)(α::Void, x::APL) = false
(==)(x::APL, α::Void) = false
(==)(α::Dict, x::APL) = false
(==)(x::APL, α::Dict) = false
(==)(α::Void, x::RationalPoly) = false
(==)(x::RationalPoly, α::Void) = false
(==)(α::Dict, x::RationalPoly) = false
(==)(x::RationalPoly, α::Dict) = false

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
eqconstant(α, t::AbstractTermLike) = α == coefficient(t) && isconstant(t)
eqconstant(t::AbstractTermLike, α) = α == coefficient(t) && isconstant(t)
eqconstant(α, p::APL) = polyeqterm(p, α)
eqconstant(p::APL, α) = polyeqterm(p, α)

function (==)(t1::AbstractTerm, t2::AbstractTerm)
    c1 = coefficient(t1)
    c2 = coefficient(t2)
    if iszero(c1)
        iszero(c2)
    else
        c1 == c2 && monomial(t1) == monomial(t2)
    end
end
(==)(p::AbstractPolynomial, t::AbstractTerm) = polyeqterm(p, t)
(==)(t::AbstractTerm, p::AbstractPolynomial) = polyeqterm(p, t)

function compare_terms(p1::AbstractPolynomial, p2::AbstractPolynomial, op)
    i1 = 1
    i2 = 1
    t1 = terms(p1)
    t2 = terms(p2)
    while true
        while i1 <= length(t1) && coefficient(t1[i1]) == 0
            i1 += 1
        end
        while i2 <= length(t2) && coefficient(t2[i2]) == 0
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
(==)(p1::AbstractPolynomial, p2::AbstractPolynomial) = compare_terms(p1, p2, ==)

(==)(p::RationalPoly, q::RationalPoly) = p.num*q.den == q.num*p.den
# Solve ambiguity with (::PolyType, ::Any)
(==)(p::APL, q::RationalPoly) = p*q.den == q.num
(==)(p, q::RationalPoly) = p*q.den == q.num

(==)(p::APL, q::MatPolynomial) = p == polynomial(q)
(==)(p::MatPolynomial, q::MatPolynomial) = iszero(p - q)

isapprox(t1::AbstractTerm, t2::AbstractTerm; kwargs...) = isapprox(coefficient(t1), coefficient(t2); kwargs...) && monomial(t1) == monomial(t2)
isapprox(p1::AbstractPolynomial, p2::AbstractPolynomial; kwargs...) = compare_terms(p1, p2, (x, y) -> isapprox(x, y; kwargs...))

function isapprox(p::MatPolynomial, q::MatPolynomial)
    p.x == q.x && isapprox(p.Q, q.Q)
end

function permcomp(f, m)
    picked = IntSet()
    for i in 1:m
        k = 0
        for j in 1:m
            if !(j in picked) && f(i, j)
                k = j
                break
            end
        end
        if k == 0
            return false
        end
        push!(picked, k)
    end
    true
end

function isapprox{C, S, T}(p::SOSDecomposition{C, S}, q::SOSDecomposition{C, T}; rtol::Real=Base.rtoldefault(S, T), atol::Real=0, ztol::Real=1e-6)
    m = length(p.ps)
    if length(q.ps) != m
        false
    else
        permcomp((i, j) -> isapprox(p.ps[i], q.ps[j], rtol=rtol, atol=atol, ztol=ztol), m)
    end
end

function isapproxzero(p::RationalPoly; ztol::Real=1e-6)
    isapproxzero(p.num, ztol=ztol)
end

isapprox(p::RationalPoly, q::RationalPoly; kwargs...) = isapprox(p.num*q.den, q.num*p.den, kwargs...)
isapprox(p::RationalPoly, q::APL; kwargs...) = isapprox(p.num, q*p.den, kwargs...)
isapprox(p::APL, q::RationalPoly; kwargs...) = isapprox(p*q.den, q.num, rtol=rtol, atol=atol, ztol=ztol)
isapprox{C}(p::RationalPoly{C}, q; kwargs...) = isapprox(p, polynomial(q), kwargs...)
isapprox{C}(p, q::RationalPoly{C}; kwargs...) = isapprox(polynomial(p), q, kwargs...)

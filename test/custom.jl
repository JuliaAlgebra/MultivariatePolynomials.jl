struct CustomTerms{T, P<:AbstractPolynomial{T}} <: AbstractPolynomialLike{T}
    p::P
end
CustomTerms(p::AbstractPolynomial{T}) where T = CustomTerms{T, typeof(p)}(p)
MultivariatePolynomials.polynomial(p::CustomTerms) = p.p

struct CustomPoly{T, P<:AbstractPolynomial{T}} <: AbstractPolynomialLike{T}
    p::P
end
CustomPoly(p::AbstractPolynomial{T}) where T = CustomPoly{T, typeof(p)}(p)
MultivariatePolynomials.terms(p::CustomPoly) = terms(p.p)

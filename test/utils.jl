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

function _typetests(x::Union{AbstractPolynomialLike{T}, Type{<:AbstractPolynomialLike{T}}}) where T
    @test (@inferred monomialtype(x)) <: AbstractMonomial

    @test (@inferred termtype(x)) <: AbstractTerm{Int}
    @test (@inferred termtype(x, Float64)) <: AbstractTerm{Float64}

    @test (@inferred polynomialtype(x)) <: AbstractPolynomial{Int}
    @test (@inferred polynomial(x)) isa AbstractPolynomial{Int}
    @test (@inferred polynomialtype(x, Float64)) <: AbstractPolynomial{Float64}
    @test (@inferred polynomial(x, Float64)) isa AbstractPolynomial{Float64}
end

function typetests(x::AbstractPolynomialLike{T}) where T
    _typetests(x)
    _typetests(typeof(x))
end

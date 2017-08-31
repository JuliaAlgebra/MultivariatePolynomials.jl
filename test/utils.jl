struct CustomPoly{T, P<:AbstractPolynomial{T}} <: AbstractPolynomialLike{T}
    p::P
end
CustomPoly(p::AbstractPolynomial{T}) where T = CustomPoly{T, typeof(p)}(p)
MultivariatePolynomials.polynomial(p::CustomPoly) = p.p

struct CustomTerms{T, P<:AbstractPolynomial{T}} <: AbstractPolynomialLike{T}
    p::P
end
CustomTerms(p::AbstractPolynomial{T}) where T = CustomTerms{T, typeof(p)}(p)
MultivariatePolynomials.terms(p::CustomTerms) = terms(p.p)

function _typetests(x, ::Type{T}) where T
    @test (@inferred coefficienttype(x)) == Int

    @test (@inferred monomialtype(x)) <: AbstractMonomial

    @test (@inferred termtype(x)) <: AbstractTerm{Int}
    @test (@inferred termtype(x, Float64)) <: AbstractTerm{Float64}

    @test (@inferred polynomialtype(x)) <: AbstractPolynomial{Int}
    @test (@inferred polynomialtype(x, Float64)) <: AbstractPolynomial{Float64}
end

function typetests(x::Union{AbstractPolynomialLike{T}, Vector{<:AbstractPolynomialLike{T}}}) where T
    _typetests(x, T)
    _typetests(typeof(x), T)
end

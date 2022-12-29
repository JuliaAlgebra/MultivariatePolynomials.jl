using Test
using MultivariatePolynomials

# Allocating size for allocating a `BigInt`.
# Half size on 32-bit.
const BIGINT_ALLOC = Sys.WORD_SIZE == 64 ? 48 : 24

function mutable_alloc_test(f, x, n)
    y = MultivariatePolynomials.MA.mutable_copy(x)
    f(y) # compile
    @test n == @allocated f(x)
end

function alloc_test(f, n)
    f() # compile
    @test n == @allocated f()
end

struct CustomPoly{T, P<:AbstractPolynomial{T}} <: AbstractPolynomialLike{T}
    p::P
end
CustomPoly(p::AbstractPolynomial{T}) where T = CustomPoly{T, typeof(p)}(p)
MultivariatePolynomials.termtype(::Type{CustomPoly{T,P}}) where {T,P} = MultivariatePolynomials.termtype(P)
MultivariatePolynomials.polynomial(p::CustomPoly) = p.p
MultivariatePolynomials.polynomial(p::CustomPoly, T::Type) = polynomial(p.p, T)
MultivariatePolynomials.variables(p::CustomPoly) = variables(p.p)
MultivariatePolynomials.monomialtype(::Type{<:CustomPoly{T, P}}) where {T, P} = monomialtype(P)
MultivariatePolynomials.constantmonomial(p::CustomPoly) = constantmonomial(p.p)
Base.copy(p::CustomPoly) = CustomPoly(copy(p.p))

struct CustomTerms{T, P<:AbstractPolynomial{T}} <: AbstractPolynomialLike{T}
    p::P
end
CustomTerms(p::AbstractPolynomial{T}) where T = CustomTerms{T, typeof(p)}(p)
MultivariatePolynomials.termtype(::Type{CustomTerms{T,P}}) where {T,P} = MultivariatePolynomials.termtype(P)
MultivariatePolynomials.terms(p::CustomTerms) = terms(p.p)
MultivariatePolynomials.variables(p::CustomTerms) = variables(p.p)
MultivariatePolynomials.monomialtype(::Type{<:CustomTerms{T, P}}) where {T, P} = monomialtype(P)
MultivariatePolynomials.constantmonomial(p::CustomPoly) = constantmonomial(p.p)
Base.copy(p::CustomTerms) = CustomTerms(copy(p.p))

function _typetests(x, ::Type{T}) where T
    @test (@inferred coefficienttype(x)) == Int

    @test (@inferred monomialtype(x)) <: AbstractMonomial

    @test (@inferred termtype(x)) <: AbstractTerm{Int}
    @test (@inferred termtype(x, Float64)) <: AbstractTerm{Float64}

    @test (@inferred polynomialtype(x)) <: AbstractPolynomial{Int}
    @test (@inferred polynomialtype(x, Float64)) <: AbstractPolynomial{Float64}

    @test (@inferred monovectype(x)) <: AbstractArray{<:AbstractMonomial}
end

function typetests(x::Union{AbstractPolynomialLike{T}, Vector{<:AbstractPolynomialLike{T}}}) where T
    _typetests(x, T)
    _typetests(typeof(x), T)
end

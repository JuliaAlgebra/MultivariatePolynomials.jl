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

struct CustomPoly{T,P<:AbstractPolynomial{T}} <: AbstractPolynomialLike{T}
    p::P
end
CustomPoly(p::AbstractPolynomial{T}) where {T} = CustomPoly{T,typeof(p)}(p)
function MultivariatePolynomials.term_type(::Type{CustomPoly{T,P}}) where {T,P}
    return MultivariatePolynomials.term_type(P)
end
MultivariatePolynomials.polynomial(p::CustomPoly) = p.p
MultivariatePolynomials.polynomial(p::CustomPoly, T::Type) = polynomial(p.p, T)
MultivariatePolynomials.variables(p::CustomPoly) = variables(p.p)
function MultivariatePolynomials.monomial_type(
    ::Type{<:CustomPoly{T,P}},
) where {T,P}
    return monomial_type(P)
end
function MultivariatePolynomials.constant_monomial(p::CustomPoly)
    return constant_monomial(p.p)
end
Base.copy(p::CustomPoly) = CustomPoly(copy(p.p))

struct CustomTerms{T,P<:AbstractPolynomial{T}} <: AbstractPolynomialLike{T}
    p::P
end
CustomTerms(p::AbstractPolynomial{T}) where {T} = CustomTerms{T,typeof(p)}(p)
function MultivariatePolynomials.term_type(::Type{CustomTerms{T,P}}) where {T,P}
    return MultivariatePolynomials.term_type(P)
end
MultivariatePolynomials.terms(p::CustomTerms) = terms(p.p)
MultivariatePolynomials.variables(p::CustomTerms) = variables(p.p)
function MultivariatePolynomials.monomial_type(
    ::Type{<:CustomTerms{T,P}},
) where {T,P}
    return monomial_type(P)
end
function MultivariatePolynomials.constant_monomial(p::CustomPoly)
    return constant_monomial(p.p)
end
Base.copy(p::CustomTerms) = CustomTerms(copy(p.p))

function _typetests(x, ::Type{T}) where {T}
    @test (@inferred coefficienttype(x)) == Int

    @test (@inferred monomial_type(x)) <: AbstractMonomial

    @test (@inferred term_type(x)) <: AbstractTerm{Int}
    @test (@inferred term_type(x, Float64)) <: AbstractTerm{Float64}

    @test (@inferred polynomial_type(x)) <: AbstractPolynomial{Int}
    @test (@inferred polynomial_type(x, Float64)) <: AbstractPolynomial{Float64}

    @test (@inferred monomial_vector_type(x)) <:
          AbstractArray{<:AbstractMonomial}
end

function typetests(
    x::Union{AbstractPolynomialLike{T},Vector{<:AbstractPolynomialLike{T}}},
) where {T}
    _typetests(x, T)
    return _typetests(typeof(x), T)
end

using Test
using LinearAlgebra
using MultivariatePolynomials
import MultivariatePolynomials as MP

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

function _typetests(x)
    @test (@inferred coefficient_type(x)) == Int

    @test (@inferred monomial_type(x))<:AbstractMonomial

    @test (@inferred term_type(x))<:AbstractTerm{Int}
    @test (@inferred term_type(x, Float64))<:AbstractTerm{Float64}

    @test (@inferred polynomial_type(x))<:AbstractPolynomial{Int}
    @test (@inferred polynomial_type(x, Float64))<:AbstractPolynomial{Float64}

    @test (@inferred monomial_vector_type(x))<:AbstractArray{<:AbstractMonomial}
end

function typetests(x)
    _typetests(x)
    return _typetests(typeof(x))
end

# A polynomial-like representation that exposes its polynomial type only.
struct CustomPolyType{T,P<:AbstractPolynomial{T}} <: AbstractPolynomialLike{T}
    p::P
end
function MultivariatePolynomials.polynomial_type(
    ::Type{CustomPolyType{T,P}},
) where {T,P}
    return P
end

# Like SOSDecomposition, promote two wrappers without calling term_type.
function Base.promote_rule(
    ::Type{CustomPolyType{T,P}},
    ::Type{CustomPolyType{S,Q}},
) where {T,P,S,Q}
    return CustomPolyType{promote_type(T, S),promote_type(P, Q)}
end

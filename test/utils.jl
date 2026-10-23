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

module MultivariatePolynomials

import LinearAlgebra

import DataStructures

import MutableArithmetics
const MA = MutableArithmetics

export AbstractPolynomialLike, AbstractTermLike, AbstractMonomialLike
"""
    AbstractPolynomialLike{T}

Abstract type for a value that can act like a polynomial. For instance, an
`AbstractTerm{T}` is an `AbstractPolynomialLike{T}` since it can act as a
polynomial of only one term.
"""
abstract type AbstractPolynomialLike{T} <: MA.AbstractMutable end

"""
    AbstractTermLike{T}

Abstract type for a value that can act like a term. For instance, an `AbstractMonomial` is an `AbstractTermLike{Int}` since it can act as a term with coefficient `1`.
"""
abstract type AbstractTermLike{T} <: AbstractPolynomialLike{T} end

"""
    AbstractMonomialLike

Abstract type for a value that can act like a monomial. For instance, an `AbstractVariable` is an `AbstractMonomialLike` since it can act as a monomial of one variable with degree `1`.
"""
abstract type AbstractMonomialLike <: AbstractTermLike{Int} end

export AbstractVariable, AbstractMonomial, AbstractTerm, AbstractPolynomial
"""
    AbstractVariable <: AbstractMonomialLike

Abstract type for a variable.
"""
abstract type AbstractVariable <: AbstractMonomialLike end

"""
    AbstractMonomial <: AbstractMonomialLike

Abstract type for a monomial, i.e. a product of variables elevated to a nonnegative integer power.
"""
abstract type AbstractMonomial <: AbstractMonomialLike end

"""
    AbstractTerm{T} <: AbstractTermLike{T}

Abstract type for a term of coefficient type `T`, i.e. the product between a value of type `T` and a monomial.
"""
abstract type AbstractTerm{T} <: AbstractTermLike{T} end

"""
    AbstractPolynomial{T} <: AbstractPolynomialLike{T}

Abstract type for a polynomial of coefficient type `T`, i.e. a sum of `AbstractTerm{T}`s.
"""
abstract type AbstractPolynomial{T} <: AbstractPolynomialLike{T} end

const APL{T} = AbstractPolynomialLike{T}

include("zip.jl")

include("variable.jl")
include("monomial.jl")
include("term.jl")
include("polynomial.jl")
include("monovec.jl")
include("ordering.jl")

include("rational.jl")

include("show.jl")
include("hash.jl")

include("promote.jl")
include("conversion.jl")

include("operators.jl")
include("comparison.jl")

include("substitution.jl")
include("differentiation.jl")
include("division.jl")
include("gcd.jl")
include("det.jl")

include("chain_rules.jl")

end # module

__precompile__()

module MultivariatePolynomials

import Base: *, +, -, /, ^, ==,
    promote_rule, convert, show, isless, size, getindex,
    one, zero, transpose, isapprox, @pure, dot, copy

abstract type AbstractPolynomialLike{T} end
abstract type AbstractTermLike{T} <: AbstractPolynomialLike{T} end
abstract type AbstractMonomialLike <: AbstractTermLike{Int} end

abstract type AbstractVariable <: AbstractMonomialLike end
abstract type AbstractMonomial <: AbstractMonomialLike end
abstract type AbstractTerm{T} <: AbstractTermLike{T} end
abstract type AbstractPolynomial{T} <: AbstractPolynomialLike{T} end

const APL{T} = AbstractPolynomialLike{T}

include("measure.jl")
include("exp.jl")

#include("utils.jl")
#include("types.jl")
#include("operators.jl")
#include("show.jl")
#include("substitution.jl")

end # module

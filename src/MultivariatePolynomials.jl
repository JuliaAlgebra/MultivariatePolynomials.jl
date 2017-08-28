__precompile__()

module MultivariatePolynomials

#using DocStringExtensions

import Base: *, +, -, /, ^, ==,
    promote_rule, convert, show, isless, size, getindex,
    one, zero, transpose, isapprox, @pure, dot, copy

export AbstractPolynomialLike, AbstractTermLike, AbstractMonomialLike
abstract type AbstractPolynomialLike{T} end
abstract type AbstractTermLike{T} <: AbstractPolynomialLike{T} end
abstract type AbstractMonomialLike <: AbstractTermLike{Int} end

export AbstractVariable, AbstractMonomial, AbstractTerm, AbstractPolynomial
abstract type AbstractVariable <: AbstractMonomialLike end
abstract type AbstractMonomial <: AbstractMonomialLike end
abstract type AbstractTerm{T} <: AbstractTermLike{T} end
abstract type AbstractPolynomial{T} <: AbstractPolynomialLike{T} end

const APL{T} = AbstractPolynomialLike{T}

include("zip.jl")
include("variable.jl")
include("monomial.jl")
include("term.jl")
include("monovec.jl")
include("polynomial.jl")

include("rational.jl")

include("conversion.jl")
include("promote.jl")
include("substitution.jl")
include("hash.jl")

include("operators.jl")
include("comparision.jl")

include("differentiation.jl")
include("division.jl")

include("show.jl")

end # module

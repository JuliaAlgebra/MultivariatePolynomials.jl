"""
    AbstractOrdering

Abstract type for a monomial ordering.
"""
abstract type AbstractOrdering end

"""
    GradedLex

Graded lexicographical monomial ordering.
Compares total monomial degree first, then breaks ties lexicographically.
Currently, a default option for all `AbstractPolynomialLike{T}`.
"""
struct GradedLex <: AbstractOrdering end

"""
    ordering(p::AbstractPolynomialLike)

Returns the ordering of polynomial `p`.
"""
ordering(::AbstractPolynomialLike) = GradedLex()

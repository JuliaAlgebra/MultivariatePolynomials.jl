"""A variable supplied by a polynomial frontend."""
abstract type AbstractVariable <: MA.AbstractMutable end

struct Variables{B,V}
    variables::V
end

abstract type AbstractMonomialIndexed end
abstract type AbstractMonomialBasis <: AbstractMonomialIndexed end

"""The monomial basis."""
struct Monomial <: AbstractMonomialBasis end

"""
    BasisElement{B,V,E}

A single element of polynomial basis `B`, identified by its variables and
exponents. `BasisElement{Monomial}` is the monomial type `AbstractMonomial`.
A coefficient times a basis element is stored in `StarAlgebras.Term`, and sums
of terms are stored in `StarAlgebras.AlgebraElement`.
"""
struct BasisElement{B<:AbstractMonomialIndexed,V,E} <: MA.AbstractMutable
    variables::Variables{B,V}
    exponents::E
end

# Noncommutative keys retain their own ordered variables and exponents. The
# parent records only the alphabet, without imposing a common factor sequence.
struct _NoncommutativeBasis{M,V} <: SA.ImplicitBasis{M,M}
    variables::V
end

const AbstractMonomial = BasisElement{Monomial}
const AbstractMonomialLike = Union{AbstractVariable,AbstractMonomial}
# Includes basis elements of non-monomial bases, e.g., Chebyshev polynomials
const _BasisElementLike = Union{AbstractVariable,BasisElement}
const _PolynomialAlgebra =
    SA.StarAlgebra{<:Union{Variables,_NoncommutativeBasis}}
const AbstractTerm{T} = SA.Term{T,<:_PolynomialAlgebra}
const AbstractPolynomial{T} = SA.AlgebraElement{T,<:_PolynomialAlgebra}
const AbstractTermLike{T} = Union{AbstractMonomialLike,AbstractTerm{T}}
const AbstractPolynomialLike{T} =
    Union{AbstractTermLike{T},AbstractPolynomial{T}}
const _APL = AbstractPolynomialLike

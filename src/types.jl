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
    Polynomial{B,V,E}

A single element of basis `B`, identified by its variables and exponents.
Sums of basis elements are stored in `StarAlgebras.AlgebraElement`.
"""
struct Polynomial{B<:AbstractMonomialIndexed,V,E} <: MA.AbstractMutable
    variables::Variables{B,V}
    exponents::E
end

# Noncommutative keys retain their own ordered variables and exponents. The
# parent records only the alphabet, without imposing a common factor sequence.
struct _NoncommutativeBasis{M,V} <: SA.ImplicitBasis{M,M}
    variables::V
end

const AbstractMonomial = Polynomial{Monomial}
const AbstractMonomialLike = Union{AbstractVariable,AbstractMonomial}
const _PolynomialAlgebra =
    SA.StarAlgebra{<:Union{Variables,_NoncommutativeBasis}}
const AbstractTerm{T} = SA.Term{T,<:_PolynomialAlgebra}
const AbstractPolynomial{T} = SA.AlgebraElement{T,<:_PolynomialAlgebra}
const AbstractTermLike{T} = Union{AbstractMonomialLike,AbstractTerm{T}}
const AbstractPolynomialLike{T} =
    Union{AbstractTermLike{T},AbstractPolynomial{T}}
const _APL = AbstractPolynomialLike

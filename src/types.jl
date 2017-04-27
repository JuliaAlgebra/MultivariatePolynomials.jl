export AbstractVariable, AbstractMonomial, AbstractTerm, AbstractPolynomial
export name, degree, vars, variables, nvars, numvariables, exponents, exponent, monomialtype, powers, coefficient, monomial, terms
export AbstractMonomialLike, AbstractTermLike, AbstractPolynomialLike

@pure vars(p) = variables(p)
@pure nvars(p) = numvariables(p)

abstract type AbstractVariable{Name} end
@pure name(::Type{<:AbstractVariable{N}}) where {N} = N
@pure name(v::AbstractVariable) = name(typeof(v))
@pure degree(::AbstractVariable) = 1
@pure numvariables(::AbstractVariable) = 1

abstract type AbstractMonomial{Variables} end
degree(m::AbstractMonomial) = sum(exponents(m))
@pure variables(::Type{<:AbstractMonomial{V}}) where {V} = V
@pure variables(m::AbstractMonomial) = variables(typeof(m))
function exponents end
function exponent end

abstract type AbstractTerm{CoeffType, MonomialType} end
monomialtype(::Type{<:AbstractTerm{C, M}}) where {C, M} = M
monomialtype(t::AbstractTerm) = monomialtype(typeof(t))
degree(t::AbstractTerm) = degree(monomial(t))
variables(T::Type{<:AbstractTerm}) = variables(monomialtype(T))
variables(t::AbstractTerm) = variables(monomialtype(t))
exponents(t::AbstractTerm) = exponents(monomial(t))
exponent(t::AbstractTerm, v::AbstractVariable) = exponent(monomial(t), v)
powers(m::AbstractMonomial) = tuplezip(variables(m), exponents(m))
function coefficient end
function monomial end

abstract type AbstractPolynomial end
# termtype(::Type{<:AbstractPolynomial{T}}) where {T} = T
# termtype(p::AbstractPolynomial) = termtype(typeof(p))
function terms end
function variables end
# variables(T::Type{<:AbstractPolynomial}) = variables(termtype(T))
# variables(p::AbstractPolynomial) = variables(termtype(p))

const AbstractMonomialLike = Union{<:AbstractVariable, <:AbstractMonomial}
const AbstractTermLike = Union{<:AbstractMonomialLike, <:AbstractTerm}
const AbstractPolynomialLike = Union{<:AbstractTermLike, <:AbstractPolynomial}

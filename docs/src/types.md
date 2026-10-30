```@meta
CurrentModule = MultivariatePolynomials
```

# Types

## Variables

```@docs
AbstractVariable
variable
name
name_base_indices
variable_union_type
similar_variable
@similar_variable
conj(::AbstractVariable)
real(::AbstractVariable)
imag(::AbstractVariable)
isreal(::AbstractVariable)
isrealpart
isimagpart
isconj
ordinary_variable
is_commutative
```

## Basis elements

A [`BasisElement`](@ref) represents one element of a chosen polynomial basis.
`BasisElement{Monomial}` is the monomial type `AbstractMonomial`.
Terms attach coefficients to basis elements, and polynomials store linear
combinations of them in `StarAlgebras.AlgebraElement`.

```@docs
BasisElement
```

## Monomials

```@docs
AbstractMonomialLike
AbstractMonomial
monomial_type
variables
effective_variables
nvariables
exponents
degree
isconstant
powers
constant_monomial
map_exponents
multiplication_preserves_monomial_order
promote_variables
promote_variables_with_maps
```

### Ordering

Use `cmp(a, b)` to compare monomials or terms with their monomial ordering.
For exponent vectors or tuples, use `cmp(order, e1, e2)`, where `order` is an
[`AbstractMonomialOrdering`](@ref) instance.

```@docs
AbstractMonomialOrdering
ordering
LexOrder
InverseLexOrder
Graded
Reverse
ExponentsIterator
```

## Terms

```@docs
AbstractTermLike
AbstractTerm
Term
term
term_type
coefficient
coefficient_type
monomial
constant_term
zero_term
degree_complex
halfdegree
```

## Polynomials

```@docs
AbstractPolynomialLike
AbstractPolynomial
polynomial
polynomial_type
terms
nterms
coefficients
coefficient(p::AbstractPolynomialLike, m::AbstractMonomialLike, vars)
monomials
mindegree
maxdegree
extdegree
leading_term
leading_coefficient
leading_monomial
deg_num_leading_terms
SA.remove_leading_term
remove_monomials
filter_terms
OfDegree
monic
map_coefficients
SA.map_coefficients!
SA.map_coefficients_to!
conj(::_APL)
real(::_APL)
imag(::_APL)
isreal(::_APL)
mindegree_complex
minhalfdegree
maxdegree_complex
maxhalfdegree
extdegree_complex
exthalfdegree
```

## Rational Polynomial Function

A rational polynomial function can be constructed with the `/` operator. Common operations such as `+`, `-`, `*`, `-` have been implemented between rational functions.
The numerator and denominator polynomials can be retrieved by the `numerator` and `denominator` functions.

## Monomial Vectors

```@docs
monomial_vector
monomial_vector_type
empty_monomial_vector
sort_monomial_vector
merge_monomial_vectors
conj(::AbstractVector{<:AbstractMonomial})
real(::AbstractVector{<:AbstractMonomial})
imag(::AbstractVector{<:AbstractMonomial})
isreal(::AbstractVector{<:AbstractMonomial})
```

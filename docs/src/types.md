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
```

### Ordering

```@docs
AbstractMonomialOrdering
ordering
compare
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
Polynomial
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
remove_leading_term
remove_monomials
filter_terms
OfDegree
monic
map_coefficients
map_coefficients!
map_coefficients_to!
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

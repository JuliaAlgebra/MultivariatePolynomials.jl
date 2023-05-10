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
ordering
mindegree
maxdegree
extdegree
leading_term
leading_coefficient
leading_monomial
remove_leading_term
remove_monomials
monic
map_coefficients
map_coefficients!
map_coefficients_to!
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
```

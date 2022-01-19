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
similarvariable
@similarvariable
```

## Monomials

```@docs
AbstractMonomialLike
AbstractMonomial
monomialtype
variables
effective_variables
nvariables
exponents
degree
isconstant
powers
constantmonomial
mapexponents
```

## Terms

```@docs
AbstractTermLike
AbstractTerm
term
termtype
coefficient
coefficienttype
monomial
constantterm
zeroterm
```

## Polynomials

```@docs
AbstractPolynomialLike
AbstractPolynomial
polynomial
polynomialtype
terms
nterms
coefficients
coefficient(p::AbstractPolynomialLike, m::AbstractMonomialLike, vars)
monomials
ordering
mindegree
maxdegree
extdegree
leadingterm
leadingcoefficient
leadingmonomial
removeleadingterm
removemonomials
monic
mapcoefficients
mapcoefficients!
mapcoefficients_to!
```

## Rational Polynomial Function

A rational polynomial function can be constructed with the `/` operator. Common operations such as `+`, `-`, `*`, `-` have been implemented between rational functions.
The numerator and denominator polynomials can be retrieved by the `numerator` and `denominator` functions.

## Monomial Vectors

```@docs
monovec
monovectype
emptymonovec
sortmonovec
mergemonovec
```

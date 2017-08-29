```@meta
CurrentModule = MultivariatePolynomials
```

# API

## Variables

```@docs
variable
name
similarvariable
@similarvariable
```

## Monomials

```@docs
monomialtype
variables
nvariables
exponents
degree
isconstant
powers
divides
constantmonomial
mapexponents
```

## Terms

```@docs
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
polynomial
polynomialtype
terms
nterms
coefficients
monomials
mindegree
maxdegree
extdegree
leadingterm
leadingcoefficient
leadingmonomial
removeleadingterm
removemonomials
monic
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

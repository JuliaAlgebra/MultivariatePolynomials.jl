# MultivariatePolynomials

[MultivariatePolynomials.jl](https://github.com/blegat/MultivariatePolynomials.jl) is an implementation independent library for manipulating multivariate polynomials.
It defines abstract types and an API for multivariate monomials, terms, polynomials, moments and measures and gives default implementation for common operations on them using the API.
If you want to manipulate multivariate polynomials easily and efficiently while being able to easily switch between different implementations, this library is exactly what you are looking for.

Supported operations are : basic arithmetic, rational polynomials, differentiation and evaluation/substitution, division and duality operations between polynomials and moments.
There is also support for solving systems of equations (soon!) and building (semi)algebraic sets.

Currently, the following implementations are available:

* [TypedPolynomials](https://github.com/rdeits/TypedPolynomials.jl)
* [DynamicPolynomials](https://github.com/blegat/DynamicPolynomials.jl)

```@contents
Pages = ["apireference.md"]
Depth = 3
```

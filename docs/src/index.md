# MultivariatePolynomials

[MultivariatePolynomials.jl](https://github.com/blegat/MultivariatePolynomials.jl) is an implementation independent library for manipulating multivariate polynomials.
It defines abstract types and an API for multivariate monomials, terms, polynomials and gives default implementation for common operations on them using the API.

On the one hand, This packages allows you to implement algorithms on multivariate polynomials that will be independant on the representation of the polynomial that will be chosen by the user.
On the other hand, it allows the user to easily switch between different representations of polynomials to see which one is faster for the algorithm that he is using.

Supported operations are : basic arithmetic, rational polynomials, evaluation/substitution, differentiation and division.

The following packages provide representations of multivariate polynomials that implement the interface:

* [TypedPolynomials](https://github.com/rdeits/TypedPolynomials.jl) : Commutative polynomials of arbitrary coefficient types
* [DynamicPolynomials](https://github.com/blegat/DynamicPolynomials.jl) : Commutative and non-commutative polynomials of arbitrary coefficient types

The following packages extend the interface and/or implement algorithms using the interface:

* [SemialgebraicSets](https://github.com/blegat/SemialgebraicSets.jl) : Sets defined by inequalities and equalities between polynomials and algorithms for solving polynomial systems of equations.
* [HomotopyContinuation](https://github.com/saschatimme/HomotopyContinuation.jl) : Solving systems of polynomials via homotopy continuation.
* [MultivariateBases](https://github.com/JuliaAlgebra/MultivariateBases.jl/) : Standardized API for multivariate polynomial bases.
* [MultivariateMoments](https://github.com/blegat/MultivariateMoments.jl) : Moments of multivariate measures and their scalar product with polynomials.
* [PolyJuMP](https://github.com/JuliaOpt/PolyJuMP.jl) : A [JuMP](https://github.com/JuliaOpt/JuMP.jl) extension for Polynomial Optimization.
* [SumOfSquares](https://github.com/JuliaOpt/SumOfSquares.jl) : Certifying the nonnegativity of polynomials, minimizing/maximizing polynomials and optimization over sum of squares polynomials using Sum of Squares Programming.

## Contents
```@contents
Pages = ["types.md", "substitution.md", "differentiation.md", "division.md", "internal.md"]
Depth = 3
```

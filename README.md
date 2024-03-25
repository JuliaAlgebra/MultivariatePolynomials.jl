# Multivariate Polynomials

| **Documentation** | **Build Status** | **Social** | **References to cite** |
|:-----------------:|:----------------:|:----------:|:----------------------:|
| [![][docs-stable-img]][docs-stable-url] | [![Build Status][build-img]][build-url] | [![Gitter][gitter-img]][gitter-url] | [![DOI][zenodo-img]][zenodo-url] |
| [![][docs-latest-img]][docs-latest-url] | [![Codecov branch][codecov-img]][codecov-url] | [<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/a/af/Discourse_logo.png/799px-Discourse_logo.png" width="64">][discourse-url] | |

**WARNING** MultivariatePolynomials v0.5 has a few breaking changes including the change in the order of monomials. This means that the leading term will now appear last in the order of monomials instead of first.

This package provides an interface for manipulating multivariate polynomials.
Implementing algorithms on polynomials using this interface will allow the algorithm to work for all polynomials implementing the interface.

The interface contains functions for accessing the coefficients, monomials, terms of the polynomial, defines arithmetic operations on them, rational functions, division with remainder, calculus/differentiation and evaluation/substitution.

## Documentation

- [**STABLE**][docs-stable-url] &mdash; **most recently tagged version of the documentation.**
- [**LATEST**][docs-latest-url] &mdash; *in-development version of the documentation.*

## Citing

Please cite the [JuliaCon 2022 presentation](https://pretalx.com/juliacon-2022/talk/TRFSJY/) [[Slides](https://drive.google.com/file/d/1q9UT5rpcmm0GdRmWm7llhVpbOx7OxGnZ/view)].
See [CITATION.bib](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl/blob/master/CITATION.bib) for the BibTeX.

## Examples

Below is a simple usage example
```julia
using TypedPolynomials
@polyvar x y # assigns x (resp. y) to a variable of name x (resp. y)
p = 2x + 3.0x*y^2 + y
@test differentiate(p, x) # compute the derivative of p with respect to x
@test differentiate.(p, (x, y)) # compute the gradient of p
@test p((x, y)=>(y, x)) # replace any x by y and y by x
@test subs(p, y=>x^2) # replace any occurence of y by x^2
@test p(x=>1, y=>2) # evaluate p at [1, 2]
```
Below is an example with `@polyvar x[1:3]`
```julia
using TypedPolynomials
A = rand(3, 3)
@polyvar x[1:3] # assign x to a tuple of variables x1, x2, x3
p = sum(x .* x) # x_1^2 + x_2^2 + x_3^2
subs(p, x[1]=>2, x[3]=>3) # x_2^2 + 13
p(x=>A*vec(x)) # corresponds to dot(A*x, A*x), need vec to convert the tuple to a vector
```

## Ecosystem

The following packages provides multivariate polynomials that implement the interface:

* [TypedPolynomials](https://github.com/rdeits/TypedPolynomials.jl) : Commutative polynomials of arbitrary coefficient types
* [DynamicPolynomials](https://github.com/JuliaAlgebra/DynamicPolynomials.jl) : Commutative and non-commutative polynomials of arbitrary coefficient types

The following packages extend/use the interface:

* [SemialgebraicSets](https://github.com/JuliaAlgebra/SemialgebraicSets.jl) : Sets defined by inequalities and equalities between polynomials and algorithms for solving polynomial systems of equations.
* [FixedPolynomials](https://github.com/JuliaAlgebra/FixedPolynomials.jl) : Fast evaluation of multivariate polynomials
* [HomotopyContinuation](https://github.com/saschatimme/HomotopyContinuation.jl) : Solving systems of polynomials via homotopy continuation.
* [MultivariateBases](https://github.com/JuliaAlgebra/MultivariateBases.jl/) : Standardized API for multivariate polynomial bases.
* [MultivariateMoments](https://github.com/JuliaAlgebra/MultivariateMoments.jl) : Moments of multivariate measures and their scalar product with polynomials.
* [PolyJuMP](https://github.com/JuliaOpt/PolyJuMP.jl) : A [JuMP](https://github.com/JuliaOpt/JuMP.jl) extension for Polynomial Optimization.
* [SumOfSquares](https://github.com/JuliaOpt/SumOfSquares.jl) : Certifying the nonnegativity of polynomials, minimizing/maximizing polynomials and optimization over sum of squares polynomials using Sum of Squares Programming.

### See also

* [Nemo](https://github.com/wbhart/Nemo.jl) for generic polynomial rings, matrix spaces, fraction fields, residue rings, power series
* [Polynomials](https://github.com/Keno/Polynomials.jl) for univariate polynomials
* [PolynomialRoots](https://github.com/giordano/PolynomialRoots.jl) for a fast complex polynomial root finder

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-stable-url]: https://JuliaAlgebra.github.io/MultivariatePolynomials.jl/stable
[docs-latest-url]: https://juliaalgebra.github.io/MultivariatePolynomials.jl/dev

[build-img]: https://github.com/JuliaAlgebra/MultivariatePolynomials.jl/workflows/CI/badge.svg?branch=master
[build-url]: https://github.com/JuliaAlgebra/MultivariatePolynomials.jl/actions?query=workflow%3ACI
[codecov-img]: http://codecov.io/github/JuliaAlgebra/MultivariatePolynomials.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/JuliaAlgebra/MultivariatePolynomials.jl?branch=master

[gitter-url]: https://gitter.im/JuliaAlgebra/Lobby?utm_source=share-link&utm_medium=link&utm_campaign=share-link
[gitter-img]: https://badges.gitter.im/JuliaAlgebra/Lobby.svg
[discourse-url]: https://discourse.julialang.org/c/domain/opt

[zenodo-url]: https://zenodo.org/badge/latestdoi/72210778
[zenodo-img]: https://zenodo.org/badge/72210778.svg

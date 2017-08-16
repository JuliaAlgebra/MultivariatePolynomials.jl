# Multivariate Polynomials

| **Documentation** | **PackageEvaluator** | **Build Status** |
|:-----------------:|:--------------------:|:----------------:|
| [![][docs-stable-img]][docs-stable-url] | [![][pkg-0.5-img]][pkg-0.5-url] | [![Build Status][build-img]][build-url] [![Build Status][winbuild-img]][winbuild-url] |
| [![][docs-latest-img]][docs-latest-url] | [![][pkg-0.6-img]][pkg-0.6-url] | [![Coveralls branch][coveralls-img]][coveralls-url] [![Codecov branch][codecov-img]][codecov-url] |

This package provides an interface for manipulating multivariate polynomials.
Implementing algorithms on polynomials using this interface will allow the algorithm to work for all polynomials implementing the interface.

The interface contains functions for accessing the coefficients, monomials, terms of the polynomial, defines arithmetic operations on them, rational functions, division with remainder, calculus/differentiation and evaluation/substitution.

The following packages provides multivariate polynomials that implement the interface.

* [TypedPolynomials](https://github.com/rdeits/TypedPolynomials.jl) : Commutative polynomials of arbitrary coefficient types
* [DynamicPolynomials](https://github.com/blegat/DynamicPolynomials.jl) : Commutative and non-commutative polynomials of arbitrary coefficient types

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
A = rand(3, 3)
@polyvar x[1:3] # assign x to a tuple of variables x1, x2, x3
p = sum(x .* x) # x_1^2 + x_2^2 + x_3^2
subs(p, x[1]=>2, x[3]=>3) # x_2^2 + 13
p(x=>A*vec(x)) # corresponds to dot(A*x, A*x), need vec to convert the tuple to a vector
```

## See also

* [Nemo](https://github.com/wbhart/Nemo.jl) for generic polynomial rings, matrix spaces, fraction fields, residue rings, power series

* [Polynomials](https://github.com/Keno/Polynomials.jl) for univariate polynomials

* [PolynomialRoots](https://github.com/giordano/PolynomialRoots.jl) for a fast complex polynomial root finder

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-stable-url]: https://blegat.github.io/MultivariatePolynomials.jl/stable
[docs-latest-url]: https://blegat.github.io/MultivariatePolynomials.jl/latest

[pkg-0.5-img]: http://pkg.julialang.org/badges/MultivariatePolynomials_0.5.svg
[pkg-0.5-url]: http://pkg.julialang.org/?pkg=MultivariatePolynomials
[pkg-0.6-img]: http://pkg.julialang.org/badges/MultivariatePolynomials_0.6.svg
[pkg-0.6-url]: http://pkg.julialang.org/?pkg=MultivariatePolynomials

[build-img]: https://travis-ci.org/blegat/MultivariatePolynomials.jl.svg?branch=master
[build-url]: https://travis-ci.org/blegat/MultivariatePolynomials.jl
[winbuild-img]: https://ci.appveyor.com/api/projects/status/4l5i8sbxev8405jl?svg=true
[winbuild-url]: https://ci.appveyor.com/project/blegat/multivariatepolynomials-jl
[coveralls-img]: https://coveralls.io/repos/github/blegat/MultivariatePolynomials.jl/badge.svg?branch=master
[coveralls-url]: https://coveralls.io/github/blegat/MultivariatePolynomials.jl?branch=master
[codecov-img]: http://codecov.io/github/blegat/MultivariatePolynomials.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/blegat/MultivariatePolynomials.jl?branch=master

# Differentiation

Given a polynomial, say ``p(x, y) = 3x^2y + x + 2y + 1``, we can differentiate it by a variable, say ``x`` and get ``\partial p(x, y) / \partial x = 6xy + 1``.
We can also differentiate it by both of its variable and get the vector ``[6xy+1, 3x^2+1]``.

```@docs
differentiate
```

# Antidifferentiation

Given a polynomial, say `p(x, y) = 3x^2y + x + 2y + 1`, we can antidifferentiate it by a variable, say `x` and get ``\int_0^x p(X, y)\mathrm{d}X = x^3y + 1/2x^2 + 2xy + x``.
We can also antidifferentiate it by both of its variable and get the vector `[x^3y + 1/2x^2 + 2xy + x, 3/2x^2y^2 + xy + y^2 + y]`.

```@docs
antidifferentiate
```

## Automatic differentiation of coefficients

Loading Mooncake enables reverse differentiation through nonmutating `+`, `-`,
`*`, scalar `/`, nonnegative integer powers, `term(c, m)`, `polynomial(c, m)`,
`polynomial(p)`, `coefficient`, and `differentiate(p, x)`. This currently covers
commutative monomial polynomials with real or complex floating-point
coefficients. Variables and monomials are symbolic constants. Other bases,
in-place operations, and direct access to coefficient storage are not covered.

```julia
using DynamicPolynomials, MultivariatePolynomials, Mooncake
@polyvar x
f(a, b, x) = coefficient(a * x + b * x, x)
cache = Mooncake.prepare_gradient_cache(f, 1.0, -1.0, x)
value, gradients = Mooncake.value_and_gradient!!(cache, f, 1.0, -1.0, x)
# value == 0.0, gradients[2:3] == (1.0, 1.0)
```

The pairing used by AD is the real part of the scalar coefficient inner
product. It is independent of the polynomial-valued `LinearAlgebra.dot`.
Polynomial cotangents have their own sparse support: cancellation of a primal
coefficient does not erase its derivative. For a polynomial argument, the
cotangent can contain monomials absent from the primal, within the argument's
variable set; it is available as `gradients[i].polynomial`. Terms have only one
coefficient coordinate. As usual with Mooncake, copy a returned gradient if it
must survive another call using the same cache.

# Antidifferentiation

Given a polynomial, say `p(x, y) = 3x^2y + x + 2y + 1`, we can antidifferentiate it by a variable, say `x` and get ``\int_0^x p(X, y)\mathrm{d}X = x^3y + 1/2x^2 + 2xy + x``.
We can also antidifferentiate it by both of its variable and get the vector `[x^3y + 1/2x^2 + 2xy + x, 3/2x^2y^2 + xy + y^2 + y]`.

```@docs
antidifferentiate
```

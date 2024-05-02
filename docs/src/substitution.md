# Subtitution

Given a polynomial, say ``p(x, y) = 3x^2y + x + 2y + 1``, one can evaluate it at a given point, e.g. ``p(2, 1) = 12 + 2 + 2 + 1 = 17`` or substitute one or more variable by a value or polynomial, e.g. ``p(x, xy^2 + 1) = 3x^2(xy^2+1) + x + 2(xy^2+1) + 1 = 3x^3y^2 + 2xy^2 + 3x^2 + x + 3``.
We distinguish the two operation as follows

* We call an evaluation an operation where **every** variable should be replace by a new value or polynomial, the syntax is `p(x => 2, y => 1)`.
* We call a subsitution an operation where **some** (or all variables) are subtituted into a new value or polynomial, the syntax is `subs(p, y => x*y^2 + 1)`.

The distinction is important for type stability for some implementations (it is important for [DynamicPolynomials](https://github.com/blegat/DynamicPolynomials.jl) but not for [TypedPolynomials](https://github.com/rdeits/TypedPolynomials.jl)).
Indeed consider a polynomial with `Int` coefficients for which we ask to replace some variables with `Int` values. If all the variables are replaced with `Int`s, the return type should be `Int`.
However, if some variables only are replaced by `Int` then the return type should be a polynomial with `Int` coefficients.

```@docs
subs
substitute
```

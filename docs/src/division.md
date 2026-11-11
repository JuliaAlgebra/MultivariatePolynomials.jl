# Division

Given two polynomials, ``p`` and ``d``, there are unique ``r`` and ``q`` such that ``p = q d + r`` and the leading term of ``d`` does not divide the leading term of ``r``.
You can obtain ``q`` using the `div` function and ``r`` using the `rem` function.
The `divrem` function returns ``(q, r)``.

Given a polynomial ``p`` and divisors ``d_1, \ldots, d_n``, one can find ``r`` and ``q_1, \ldots, q_n`` such that ``p = q_1 d_1 + \cdots + q_n d_n + r`` and none of the leading terms of ``q_1, \ldots, q_n`` divide the leading term of ``r``.
You can obtain the vector ``[q_1, \ldots, q_n]`` using `div(p, d)` where ``d = [d_1, \ldots, d_n]`` and ``r`` using the `rem` function with the same arguments.
The `divrem` function returns ``(q, r)``.
```@docs
divides
div_multiple
```

Note that the coefficients of the polynomials need to be a field for `div`,
`rem` and `divrem` to work.
If the coefficient type is not a field, it is promoted to a field using [`promote_to_field`](@ref).
```@docs
promote_to_field
```

For noncommutative polynomials, `rem(p, d)` and `rem(p, [d_1, ..., d_n])`
reduce contiguous occurrences of each divisor's leading monomial, preserving
the factors on both sides. For example, with noncommutative variables:

```julia
rem(x * y^2 * z, y^2 - (2x - 1)) == 2x^2 * z - x * z
```

Divisors are tried in the supplied order, and the leftmost occurrence is
reduced first. Termination requires a well-founded monomial order compatible
with multiplication on both sides. The remainder can depend on these choices
unless the divisors form a noncommutative Gröbner basis for that order;
`rem` does not compute or complete such a basis.
This two-sided reduction is supported only by `rem`, since its reduction
steps cannot generally be represented by the left quotients returned by
`div` and `divrem`.

Alternatively, [`pseudo_rem`](@ref) or [`pseudo_divrem`](@ref) can be used
instead as they do not require the coefficient type to be a field.
```@docs
pseudo_rem
pseudo_divrem
rem_or_pseudo_rem
```

## Greatest Common Divisor (GCD)

The Greatest Common Divisor (GCD) and Least Common Multiple (LCM) can be
obtained for integers respectively with the `gcd` and `lcm` functions.
The same functions can be used with monomials and polynomials:
```@docs
gcd
lcm
AbstractUnivariateGCDAlgorithm
GeneralizedEuclideanAlgorithm
SubresultantAlgorithm
```
Internal functions of the `gcd` algorithm:
```@docs
isolate_variable
primitive_univariate_gcd!
univariate_gcd
content
primitive_part
primitive_part_content
```

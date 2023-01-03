# Division

The `gcd` and `lcm` functions of `Base` have been implemented for monomials, you have for example `gcd(x^2*y^7*z^3, x^4*y^5*z^2)` returning `x^2*y^5*z^2` and `lcm(x^2*y^7*z^3, x^4*y^5*z^2)` returning `x^4*y^7*z^3`.

Given two polynomials, ``p`` and ``d``, there are unique ``r`` and ``q`` such that ``p = q d + r`` and the leading term of ``d`` does not divide the leading term of ``r``.
You can obtain ``q`` using the `div` function and ``r`` using the `rem` function.
The `divrem` function returns ``(q, r)``.

Given a polynomial ``p`` and divisors ``d_1, \ldots, d_n``, one can find ``r`` and ``q_1, \ldots, q_n`` such that ``p = q_1 d_1 + \cdots + q_n d_n + r`` and none of the leading terms of ``q_1, \ldots, q_n`` divide the leading term of ``r``.
You can obtain the vector ``[q_1, \ldots, q_n]`` using `div(p, d)` where ``d = [d_1, \ldots, d_n]`` and ``r`` using the `rem` function with the same arguments.
The `divrem` function returns ``(q, r)``.

```@docs
divides
pseudo_rem
rem_or_pseudo_rem
```

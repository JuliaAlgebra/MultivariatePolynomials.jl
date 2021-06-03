export variables, nvariables, exponents, degree, isconstant, powers, constantmonomial, mapexponents

"""
    monomialtype(p::AbstractPolynomialLike)

Return the type of the monomials of `p`.

    termtype(::Type{PT}) where PT<:AbstractPolynomialLike

Returns the type of the monomials of a polynomial of type `PT`.
"""
monomialtype(::Union{M, Type{M}}) where M<:AbstractMonomial = M
monomialtype(::Union{PT, Type{PT}}) where PT <: APL = monomialtype(termtype(PT))
monomialtype(::Union{AbstractVector{PT}, Type{<:AbstractVector{PT}}}) where PT <: APL = monomialtype(PT)

"""
    variables(p::AbstractPolynomialLike)

Returns the tuple of the variables of `p` in decreasing order. It could contain variables of zero degree, see the example section.

### Examples

Calling `variables(x^2*y)` should return `(x, y)` and calling `variables(x)` should return `(x,)`.
Note that the variables of `m` does not necessarily have nonzero degree.
For instance, `variables([x^2*y, y*z][1])` is usually `(x, y, z)` since the two monomials have been promoted to a common type.
"""
function variables end
variables(t::AbstractTerm) = variables(monomial(t))

"""
    nvariables(p::AbstractPolynomialLike)

Returns the number of variables in `p`, i.e. `length(variables(p))`. It could be more than the number of variables with nonzero degree (see the Examples section of [`variables`](@ref)).

### Examples

Calling `nvariables(x^2*y)` should return at least 2 and calling `nvariables(x)` should return at least 1.
"""
nvariables(::Union{AbstractVariable, Type{<:AbstractVariable}}) = 1
nvariables(t::AbstractTerm) = nvariables(monomial(t))
nvariables(p::APL) = length(variables(p))

"""
    exponents(t::AbstractTermLike)

Returns the exponent of the variables in the monomial of the term `t`.

### Examples

Calling `exponents(x^2*y)` should return `(2, 1)`.
"""
exponents(t::AbstractTerm) = exponents(monomial(t))
exponents(v::AbstractVariable) = (1,)

"""
    degree(t::AbstractTermLike)

Returns the *total degree* of the monomial of the term `t`, i.e. `sum(exponents(t))`.

    degree(t::AbstractTermLike, v::AbstractVariable)

Returns the exponent of the variable `v` in the monomial of the term `t`.

### Examples

Calling `degree(x^2*y)` should return 3 which is ``2 + 1``.
Calling `degree(x^2*y, x)` should return 2 and calling `degree(x^2*y, y)` should return 1.

"""
degree(t::AbstractTermLike) = sum(exponents(t))

degree(t::AbstractTermLike, var::AbstractVariable) = degree(monomial(t), var)
degree(v::AbstractVariable, var::AbstractVariable) = (v == var ? 1 : 0)
#_deg(v::AbstractVariable) = 0
#_deg(v::AbstractVariable, power, powers...) = v == power[1] ? power[2] : _deg(v, powers...)
#degree(m::AbstractMonomial, v::AbstractVariable) = _deg(v, powers(t)...)

function degree(m::AbstractMonomial, v::AbstractVariable)
    deg = 0
    # With Noncommutative variables, there may be several powers with the same variable
    for (var, exp) in powers(m)
        if var == v
            deg += exp
        end
    end
    return deg
end

"""
    isconstant(t::AbstractTermLike)

Returns whether the monomial of `t` is constant.
"""
isconstant(t::AbstractTermLike) = all(iszero.(exponents(t)))
isconstant(v::AbstractVariable) = false

"""
    powers(t::AbstractTermLike)

Returns an iterator over the powers of the monomial of `t`.

### Examples

Calling `powers(3x^4*y) should return `((x, 4), (y, 1))`.
"""
powers(t::AbstractTermLike) = _zip(variables(t), exponents(t))

"""
    constantmonomial(p::AbstractPolynomialLike)

Returns a constant monomial of the monomial type of `p` with the same variables as `p`.

    constantmonomial(::Type{PT}) where {PT<:AbstractPolynomialLike}

Returns a constant monomial of the monomial type of a polynomial of type `PT`.
"""
function constantmonomial end

"""
    mapexponents(f, m1::AbstractMonomialLike, m2::AbstractMonomialLike)

If ``m_1 = \\prod x^{\\alpha_i}`` and ``m_2 = \\prod x^{\\beta_i}`` then it returns the monomial ``m = \\prod x^{f(\\alpha_i, \\beta_i)}``.

### Examples

The multiplication `m1 * m2` is equivalent to `mapexponents(+, m1, m2)`, the unsafe division `_div(m1, m2)` is equivalent to `mapexponents(-, m1, m2)`, `gcd(m1, m2)` is equivalent to `mapexponents(min, m1, m2)`, `lcm(m1, m2)` is equivalent to `mapexponents(max, m1, m2)`.
"""
mapexponents(f, m1::AbstractMonomialLike, m2::AbstractMonomialLike) = mapexponents(f, monomial(m1), monomial(m2))

function mapexponents_to! end
function mapexponents! end

Base.one(::Type{TT}) where {TT<:AbstractMonomialLike} = constantmonomial(TT)
Base.one(t::AbstractMonomialLike) = constantmonomial(t)
MA.promote_operation(::typeof(one), MT::Type{<:AbstractMonomialLike}) = monomialtype(MT)
# See https://github.com/JuliaAlgebra/MultivariatePolynomials.jl/issues/82
# By default, Base do oneunit(v::VT) = VT(one(v)).
# This tries to convert a monomial to a variable which does not work.
# The issue here is there is no way to represent the multiplicative identity
# using the variable type. The best we can do is return a monomial even
# if it does not exactly match the definition of oneunit.
Base.oneunit(v::AbstractVariable) = one(v)
Base.oneunit(VT::Type{<: AbstractVariable}) = one(VT)

"""
    monomial_type(p::AbstractPolynomialLike)

Return the type of the monomials of `p`.

    term_type(::Type{PT}) where PT<:AbstractPolynomialLike

Returns the type of the monomials of a polynomial of type `PT`.
"""
monomial_type(::Union{M,Type{M}}) where {M<:AbstractMonomial} = M
function monomial_type(::Union{PT,Type{PT}}) where {PT<:_APL}
    return monomial_type(term_type(PT))
end
function monomial_type(
    ::Union{AbstractVector{PT},Type{<:AbstractVector{PT}}},
) where {PT<:_APL}
    return monomial_type(PT)
end

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
function variables(
    ::Type{PT},
) where {PT<:Union{AbstractPolynomial,AbstractTerm}}
    return variables(monomial_type(PT))
end

"""
    nvariables(p::AbstractPolynomialLike)

Returns the number of variables in `p`, i.e. `length(variables(p))`. It could be more than the number of variables with nonzero degree (see the Examples section of [`variables`](@ref)).

### Examples

Calling `nvariables(x^2*y)` should return at least 2 and calling `nvariables(x)` should return at least 1.
"""
nvariables(::Union{AbstractVariable,Type{<:AbstractVariable}}) = 1
nvariables(t::AbstractTerm) = nvariables(monomial(t))
nvariables(::Type{TT}) where {TT<:AbstractTerm} = variables(monomial_type(TT))
nvariables(p::_APL) = length(variables(p))

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
degree(t::AbstractTermLike) = sum(exponents(t); init = 0)

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
isconstant(t::AbstractTermLike) = all(iszero, exponents(t))
isconstant(v::AbstractVariable) = false

"""
    powers(t::AbstractTermLike)

Returns an iterator over the powers of the monomial of `t`.

### Examples

Calling `powers(3x^4*y) should return `((x, 4), (y, 1))`.
"""
powers(t::AbstractTermLike) = _zip(variables(t), exponents(t))

"""
    constant_monomial(p::AbstractPolynomialLike)

Returns a constant monomial of the monomial type of `p` with the same variables as `p`.

    constant_monomial(::Type{PT}) where {PT<:AbstractPolynomialLike}

Returns a constant monomial of the monomial type of a polynomial of type `PT`.
"""
function constant_monomial end
function constant_monomial(::Type{MT}) where {MT<:AbstractMonomial}
    return error("`constant_monomial` not implemented for $MT.")
end
function constant_monomial(::Type{PT}) where {PT<:_APL}
    return constant_monomial(monomial_type(PT))
end
constant_monomial(t::AbstractTerm) = constant_monomial(monomial(t))

"""
    map_exponents(f, m1::AbstractMonomialLike, m2::AbstractMonomialLike)

If ``m_1 = \\prod x^{\\alpha_i}`` and ``m_2 = \\prod x^{\\beta_i}`` then it returns the monomial ``m = \\prod x^{f(\\alpha_i, \\beta_i)}``.

### Examples

The multiplication `m1 * m2` is equivalent to `map_exponents(+, m1, m2)`, the unsafe division `div_multiple(m1, m2)` is equivalent to `map_exponents(-, m1, m2)`, `gcd(m1, m2)` is equivalent to `map_exponents(min, m1, m2)`, `lcm(m1, m2)` is equivalent to `map_exponents(max, m1, m2)`.
"""
function map_exponents(f, m1::AbstractMonomialLike, m2::AbstractMonomialLike)
    return map_exponents(f, monomial(m1), monomial(m2))
end

function map_exponents_to! end
function map_exponents! end

map_exponents(f, a, b, ::MA.IsMutable) = map_exponents!(f, a, b)
map_exponents(f, a, b, ::MA.IsNotMutable) = map_exponents(f, a, b)

Base.one(::Type{TT}) where {TT<:AbstractMonomialLike} = constant_monomial(TT)
Base.one(t::AbstractMonomialLike) = constant_monomial(t)
function MA.promote_operation(::typeof(one), MT::Type{<:AbstractMonomialLike})
    return monomial_type(MT)
end
# See https://github.com/JuliaAlgebra/MultivariatePolynomials.jl/issues/82
# By default, Base do oneunit(v::VT) = VT(one(v)).
# This tries to convert a monomial to a variable which does not work.
# The issue here is there is no way to represent the multiplicative identity
# using the variable type. The best we can do is return a monomial even
# if it does not exactly match the definition of oneunit.
Base.oneunit(v::AbstractVariable) = one(v)
Base.oneunit(VT::Type{<:AbstractVariable}) = one(VT)

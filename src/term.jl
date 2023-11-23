"""
    term(coef, mono::AbstractMonomialLike)

Returns a term with coefficient `coef` and monomial `mono`. There are two key
difference between this and `coef * mono`:

* `term(coef, mono)` does not copy `coef` and `mono` so modifying this term
  with MutableArithmetics may modifying the input of this function.
  To avoid this, call `term(MA.copy_if_mutable(coef), MA.copy_if_mutable(mono))`
  where `MA = MutableArithmetics`.
* Suppose that `coef = (x + 1)` and `mono = x^2`, `coef * mono` gives the polynomial
  with integer coefficients `x^3 + x^2` which `term(x + 1, x^2)` gives a term
  with polynomial coefficient `x + 1`.

    term(p::AbstractPolynomialLike)

Converts the polynomial `p` to a term.
When applied on a polynomial, it throws an `InexactError` if it has more than one term.
When applied to a term, it is the identity and does not copy it.
When applied to a monomial, it create a term of type `AbstractTerm{Bool}`.
"""
function term end
term(coef, var::AbstractVariable) = term(coef, monomial(var))
function term(coef, mono::AbstractMonomialLike)
    return term_type(mono, typeof(coef))(coef, mono)
end
term(p::_APL) = convert(term_type(typeof(p)), p)

"""
    term_type(p::AbstractPolynomialLike)

Returns the type of the terms of `p`.

    term_type(::Type{PT}) where PT<:AbstractPolynomialLike

Returns the type of the terms of a polynomial of type `PT`.

    term_type(p::AbstractPolynomialLike, ::Type{T}) where T

Returns the type of the terms of `p` but with coefficient type `T`.

    term_type(::Type{PT}, ::Type{T}) where {PT<:AbstractPolynomialLike, T}

Returns the type of the terms of a polynomial of type `PT` but with coefficient type `T`.
"""
term_type(::Type{T}) where {T<:AbstractTerm} = T
term_type(p::Type{<:_APL}, ::Type{T}) where {T} = term_type(term_type(p), T)
term_type(::Type{M}) where {M<:AbstractMonomialLike} = term_type(M, Bool)
term_type(v::Type{<:AbstractVariable}) = term_type(monomial_type(v))
function term_type(v::Type{<:AbstractVariable}, ::Type{T}) where {T}
    return term_type(monomial_type(v), T)
end
term_type(p::_APL, ::Type{T}) where {T} = term_type(typeof(p), T)
term_type(p::_APL) = term_type(typeof(p))
function term_type(
    ::Union{AbstractVector{PT},Type{<:AbstractVector{PT}}},
) where {PT<:_APL}
    return term_type(PT)
end
function term_type(
    ::Union{AbstractVector{PT},Type{<:AbstractVector{PT}}},
    ::Type{T},
) where {PT<:_APL,T}
    return term_type(PT, T)
end

"""
    coefficient(t::AbstractTermLike)

Returns the coefficient of the term `t`.

    coefficient(p::AbstractPolynomialLike, m::AbstractMonomialLike)

Returns the coefficient of the monomial `m` in `p`.

### Examples

Calling `coefficient` on ``4x^2y`` should return ``4``.
Calling `coefficient(2x + 4y^2 + 3, y^2)` should return ``4``.
Calling `coefficient(2x + 4y^2 + 3, x^2)` should return ``0``.
"""
function coefficient end
coefficient(t::AbstractTerm) = t.coefficient # by convention, the field should be `coefficient`
coefficient(m::AbstractMonomialLike) = 1
function coefficient(
    p::AbstractPolynomialLike{T},
    m::AbstractMonomialLike,
) where {T}
    for t in terms(p)
        if monomial(t) == m
            return coefficient(t)
        end
    end
    return zero(T)
end

"""
    coefficient(p::AbstractPolynomialLike, m::AbstractMonomialLike, vars)::AbstractPolynomialLike

Returns the coefficient of the monomial `m` of the polynomial `p` considered as a polynomial in variables
`vars`.

### Example
Calling `coefficient((a+b)x^2+2x+y*x^2, x^2, [x,y])` should return `a+b`.
Calling `coefficient((a+b)x^2+2x+y*x^2, x^2, [x])` should return `a+b+y`.
"""
function coefficient(f::_APL, m::AbstractMonomialLike, vars)
    coeff = zero(f)
    for t in terms(f)
        match = true
        for v in vars
            if degree(t, v) != degree(m, v)
                match = false
                break
            end
        end
        match || continue

        coeff += term(coefficient(t), div_multiple(monomial(t), m))
    end
    return coeff
end

"""
    coefficient_type(p::AbstractPolynomialLike)

Returns the coefficient type of `p`.

    coefficient_type(::Type{PT}) where PT

Returns the coefficient type of a polynomial of type `PT`.

### Examples

Calling `coefficient_type` on ``(4//5)x^2y`` should return `Rational{Int}`,
calling `coefficient_type` on ``1.0x^2y + 2.0x`` should return `Float64` and
calling `coefficient_type` on ``xy`` should return `Bool`.
"""
function coefficient_type(
    ::Union{PT,Type{PT},AbstractVector{PT},Type{<:AbstractVector{PT}}},
) where {T,PT<:_APL{T}}
    return T
end
#coefficient_type(::{T, Type{T}}) where {T} = T

"""
    monomial(t::AbstractTermLike)

Returns the monomial of the term `t`.

### Examples

Calling `monomial` on ``4x^2y`` should return ``x^2y``.
"""
function monomial end
monomial(t::AbstractTerm) = t.monomial # by convention, the field should be `monomial`.
monomial(m::AbstractMonomial) = m

"""
    constant_term(α, p::AbstractPolynomialLike)

Creates a constant term with coefficient α and the same variables as p.

    constant_term(α, ::Type{PT} where {PT<:AbstractPolynomialLike}

Creates a constant term of the term type of a polynomial of type `PT`.
"""
constant_term(α, p) = term(α, constant_monomial(p))

# zero should return a polynomial since it is often used to keep the result of a summation of terms.
# For example, Base.vecdot(x::Vector{<:AbstractTerm}, y:Vector{Int}) starts with `s = zero(dot(first(x), first(y)))` and then adds terms.
# We want `s` to start as a polynomial for this operation to be efficient.
#Base.zero(::Type{TT}) where {T, TT<:AbstractTermLike{T}} = zero(T) * constant_monomial(TT)
#Base.zero(t::AbstractTermLike{T}) where {T} = zero(T) * constant_monomial(t)
"""
    zero_term(p::AbstractPolynomialLike{T}) where T

Equivalent to `constant_term(zero(T), p)`.

    zero_term(α, ::Type{PT} where {T, PT<:AbstractPolynomialLike{T}}

Equivalent to `constant_term(zero(T), PT)`.
"""
zero_term(::Type{PT}) where {T,PT<:_APL{T}} = constant_term(zero(T), PT)
zero_term(p::_APL{T}) where {T} = constant_term(zero(T), p)

function Base.zero(::Type{TT}) where {T,TT<:AbstractTermLike{T}}
    return zero(polynomial_type(TT))
end
Base.zero(t::AbstractTermLike{T}) where {T} = zero(polynomial_type(t))
function MA.promote_operation(::typeof(zero), PT::Type{<:AbstractTermLike})
    return polynomial_type(PT)
end
function Base.one(::Type{TT}) where {T,TT<:AbstractTermLike{T}}
    return term(one(T), constant_monomial(TT))
end
Base.one(t::AbstractTermLike{T}) where {T} = term(one(T), constant_monomial(t))
function MA.promote_operation(::typeof(one), TT::Type{<:AbstractTermLike})
    return term_type(TT)
end
function MA.promote_operation(::typeof(one), PT::Type{<:AbstractPolynomialLike})
    return polynomial_type(PT)
end

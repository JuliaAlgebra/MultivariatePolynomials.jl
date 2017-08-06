export polynomial, polynomialtype, terms, nterms, coefficients, monomials
export coefficienttype, monomialtype
export mindeg, maxdeg, extdeg
export leadingterm, leadingcoefficient, leadingmonomial
export removeleadingterm, removemonomials
export vars, nvars

Base.norm(p::AbstractPolynomialLike, r::Int=2) = norm(coefficients(p), r)

function Base.hash(p::AbstractPolynomial, u::UInt)
    if iszero(p)
        hash(0, u)
    else
        reduce((u, t) -> hash(t, u), u, terms(p))
    end
end

# Conversion polynomial -> scalar
function Base.convert{S}(::Type{S}, p::APL)
    s = zero(S)
    for t in p
        if !isconstant(t)
            # The polynomial is not constant
            throw(InexactError())
        end
        s += S(coefficient(t))
    end
    s
end
Base.convert{PT<:APL}(::Type{PT}, p::PT) = p

coefficienttype{T}(::Type{<:APL{T}}) = T
coefficienttype{T}(::APL{T}) = T
#coefficienttype{T}(::Type{T}) = T
#coefficienttype{T}(::T) = T

monomialtype{T}(::Type{<:APL{T}}) = T
monomialtype{T}(::APL{T}) = T

changecoefficienttype{TT<:AbstractTermLike, T}(::Type{TT}, ::Type{T}) = termtype(TT, T)
changecoefficienttype{PT<:AbstractPolynomial, T}(::Type{PT}, ::Type{T}) = polynomialtype(PT, T)

changecoefficienttype{PT<:APL, T}(p::PT, ::Type{T}) = changecoefficienttype(PT, T)(p)

"""
    polynomial(p::AbstractPolynomialLike)

Converts `p` to a value with polynomial type.

    polynomial(a::AbstractVector, mv::AbstractVector{<:AbstractMonomialLike})

Creates a polynomial equal to `dot(a, mv)`.

    polynomial(terms::AbstractVector{<:AbstractTerm})

Creates a polynomial equal to `sum(terms)`.

### Examples

Calling `polynomial([2, 4, 1], [x, x^2*y, x*y])` should return ``4x^2y + xy + 2x``.
"""
polynomial(ts::AbstractVector{<:AbstractTerm}) = polynomial(coefficient.(ts), monomial.(ts))

polynomial{T}(Q::AbstractMatrix{T}, mv::AbstractVector) = polynomial(Q, mv, T)

polynomialtype(p::APL) = Base.promote_op(polynomial, p)

"""
    terms(p::AbstractPolynomialLike)

Returns an iterator over the nonzero terms of the polynomial `p` sorted in the decreasing monomial order.

### Examples

Calling `terms` on ``4x^2y + xy + 2x`` should return an iterator of ``[4x^2y, xy, 2x]``.
"""
terms(t::AbstractTermLike) = (term(t),)

"""
    nterms(p::AbstractPolynomialLike)

Returns the number of nonzero terms in `p`, i.e. `length(terms(p))`.

### Examples

Calling `nterms` on ``4x^2y + xy + 2x`` should return 3.
"""
function nterms end

nterms(::AbstractTermLike) = 1
nterms(p::AbstractPolynomialLike) = length(terms(p))

"""
    coefficients(p::AbstractPolynomialLike)

Returns an iterator over the coefficients of `p` of the nonzero terms of the polynomial sorted in the decreasing monomial order.

### Examples

Calling `coefficients` on ``4x^2y + xy + 2x`` should return an iterator of ``[4, 1, 2]``.
"""
function coefficients end

"""
    monomials(p::AbstractPolynomialLike)

Returns an iterator over the monomials of `p` of the nonzero terms of the polynomial sorted in the decreasing order.

    monomials(vars::Tuple, degs::AbstractVector{Int}, filter::Function = m -> true)

Builds the vector of all the monovec `m` with variables `vars` such that the degree `deg(m)` is in `degs` and `filter(m)` is `true`.

### Examples

Calling `monomials` on ``4x^2y + xy + 2x`` should return an iterator of ``[x^2y, xy, x]``.

Calling `monomials((x, y), [1, 3], m -> exponent(m, y) != 1)` should return `[x^3, x*y^2, y^3, x]` where `x^2*y` and `y` have been excluded by the filter.
"""
function monomials end

"""
$(SIGNATURES)

Returns the minimal total degree of the monomials of `p`, i.e. `minimum(deg, terms(p))`.

### Examples
Calling `mindeg` on on ``4x^2y + xy + 2x`` should return 1.
"""
function mindeg(p::AbstractPolynomialLike)
    minimum(deg, terms(p))
end

"""
$(SIGNATURES)

Returns the maximal total degree of the monomials of `p`, i.e. `maximum(deg, terms(p))`.

### Examples
Calling `maxdeg` on on ``4x^2y + xy + 2x`` should return 3.
"""
function maxdeg(p::AbstractPolynomialLike)
    maximum(deg, terms(p))
end

"""
$(SIGNATURES)

Returns the extremal total degrees of the monomials of `p`, i.e. `(mindeg(p), maxdeg(p))`.

### Examples
Calling `extdeg` on on ``4x^2y + xy + 2x`` should return `(1, 3)`.
"""
function extdeg(p::AbstractPolynomialLike)
    (mindeg(p), maxdeg(p))
end

"""
    leadingterm(p::AbstractPolynomialLike)

Returns the coefficient of the leading term, i.e. `first(terms(p))`.

### Examples

Calling `leadingterm` on ``4x^2y + xy + 2x`` should return ``4x^2y``.
"""
function leadingterm(p::AbstractPolynomialLike)
    first(terms(p))
end
leadingterm(t::AbstractTermLike) = term(t)

"""
$(SIGNATURES)

Returns the coefficient of the leading term of `p`, i.e. `coefficient(leadingterm(p))`.

### Examples

Calling `leadingcoefficient` on ``4x^2y + xy + 2x`` should return ``4`` and calling it on ``0`` should return ``0``.
"""
function leadingcoefficient(p::AbstractPolynomialLike)
    coefficient(leadingterm(p))
end

"""
$(SIGNATURES)

Returns the monomial of the leading term of `p`, i.e. `monomial(leadingterm(p))` or `first(monomials(p))`.

### Examples

Calling `leadingmonomial` on ``4x^2y + xy + 2x`` should return ``x^2y``.
"""
function leadingmonomial(p::AbstractPolynomialLike)
    # first(monomials(p)) would be more efficient for DynamicPolynomials but
    # monomial(leadingterm(p)) is more efficient for TypedPolynomials and is better if p is a term
    monomial(leadingterm(p))
end

"""
$(SIGNATURES)

Returns a polynomial with the leading term removed in the polynomial `p`.

### Examples

Calling `removeleadingterm` on ``4x^2y + xy + 2x`` should return ``xy + 2x``.
"""
function removeleadingterm(p::AbstractPolynomialLike)
    polynomial(Iterators.drop(terms(p), 1))
end

#$(SIGNATURES)
"""

Returns a polynomial with the terms having their monomial in the monomial vector `mv` removed in the polynomial `p`.

### Examples

Calling `removemonomials(4x^2*y + x*y + 2x, [x*y])` should return ``4x^2*y + 2x``.
"""
function removemonomials{MT <: AbstractMonomialLike}(p::AbstractPolynomialLike, mv::AbstractVector{MT})
    smv = monovec(mv) # Make sure it is sorted
    i = 1
    q = zero(p)
    for t in terms(p)
        m = monomial(t)
        while i <= length(smv) && smv[i] > m
            i += 1
        end
        if i <= length(smv) && smv[i] == m
            q += t
        end
    end
    q
end

"""
    vars(p::AbstractPolynomialLike)

Returns the tuple of the variables of `p` in decreasing order. It could contain variables of zero degree, see the example section.

### Examples

Calling `vars(x^2*y)` should return `(x, y)` and calling `vars(x)` should return `(x,)`.
Note that the variables of `m` does not necessarily have nonzero exponent.
For instance, `vars([x^2*y, y*z][1])` is usually `(x, y, z)` since the two monomials have been promoted to a common type.
"""
function vars end

vars(x::AbstractVariable) = (x,)

"""
    nvars(p::AbstractPolynomialLike)

Returns the number of variables in `p`, i.e. `length(vars(p))`. It could be more than the number of variables with nonzero exponent (see [the Examples section of `vars`](@ref MultivariatePolynomials.vars)).

### Examples

Calling `nvars(x^2*y)` should return at least 2 and calling `nvars(x)` should return at least 1.
"""
function nvars(p::AbstractPolynomialLike)
    length(vars(p))
end

nvars(::AbstractVariable) = 1

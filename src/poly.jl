export polynomial, polynomialtype, terms, nterms, coefficients, monomials
export coefficienttype, monomialtype
export mindeg, maxdeg, extdeg
export leadingterm, leadingcoefficient, leadingmonomial
export removeleadingterm, removemonomials
export variables, nvariables

Base.norm(p::AbstractPolynomialLike, r::Int=2) = norm(coefficients(p), r)

function Base.hash(p::AbstractPolynomial, u::UInt)
    if iszero(p)
        hash(0, u)
    else
        reduce((u, t) -> hash(t, u), u, terms(p))
    end
end

Base.convert(::Type{Any}, p::APL) = p
# Conversion polynomial -> scalar
function Base.convert(::Type{S}, p::APL) where {S}
    s = zero(S)
    for t in terms(p)
        if !isconstant(t)
            # The polynomial is not constant
            throw(InexactError())
        end
        s += S(coefficient(t))
    end
    s
end
Base.convert(::Type{PT}, p::PT) where {PT<:APL} = p
function Base.convert(::Type{MT}, t::AbstractTerm) where {MT<:AbstractMonomial}
    if isone(coefficient(t))
        monomial(t)
    else
        error("Cannot convert a term with a coefficient that is not one into a monomial")
    end
end

coefficienttype(::Type{<:APL{T}}) where {T} = T
coefficienttype(::APL{T}) where {T} = T
#coefficienttype(::Type{T}) where {T} = T
#coefficienttype(::T) where {T} = T

monomialtype(::Type{M}) where M<:AbstractMonomial = M

changecoefficienttype(::Type{TT}, ::Type{T}) where {TT<:AbstractTermLike, T} = termtype(TT, T)
changecoefficienttype(::Type{PT}, ::Type{T}) where {PT<:AbstractPolynomial, T} = polynomialtype(PT, T)

changecoefficienttype(p::PT, ::Type{T}) where {PT<:APL, T} = changecoefficienttype(PT, T)(p)

"""
    polynomial(p::AbstractPolynomialLike)

Converts `p` to a value with polynomial type.

    polynomial(p::AbstractPolynomialLike, ::Type{T}) where T

Converts `p` to a value with polynomial type with coefficient type `T`.

    polynomial(a::AbstractVector, mv::AbstractVector{<:AbstractMonomialLike})

Creates a polynomial equal to `dot(a, mv)`.

    polynomial(terms::AbstractVector{<:AbstractTerm})

Creates a polynomial equal to `sum(terms)`.

    polynomial(f::Function, mv::AbstractVector{<:AbstractMonomialLike})

Creates a polynomial equal to `sum(f(i) * mv[i] for i in 1:length(mv))`.

### Examples

Calling `polynomial([2, 4, 1], [x, x^2*y, x*y])` should return ``4x^2y + xy + 2x``.
"""
polynomial(p::AbstractPolynomial) = p
polynomial(ts::AbstractVector{<:AbstractTerm}) = polynomial(coefficient.(ts), monomial.(ts))
polynomial(a::AbstractVector, x::AbstractVector) = polynomial([α * m for (α, m) in zip(a, x)])
polynomial(f::Function, mv::AbstractVector{<:AbstractMonomialLike}) = polynomial([f(i) * mv[i] for i in 1:length(mv)])
function polynomial(Q::AbstractMatrix, mv::AbstractVector)
    dot(mv, Q * mv)
end
function polynomial(Q::AbstractMatrix, mv::AbstractVector, ::Type{T}) where T
    polynomial(polynomial(Q, mv), T)
end

polynomialtype(p::APL) = Base.promote_op(polynomial, p)

"""
    terms(p::AbstractPolynomialLike)

Returns an iterator over the nonzero terms of the polynomial `p` sorted in the decreasing monomial order.

### Examples

Calling `terms` on ``4x^2y + xy + 2x`` should return an iterator of ``[4x^2y, xy, 2x]``.
"""
terms(t::AbstractTermLike) = (term(t),)
terms(p::AbstractPolynomialLike) = terms(polynomial(p))

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
coefficients(p::APL) = coefficient.(terms(p))

"""
    monomials(p::AbstractPolynomialLike)

Returns an iterator over the monomials of `p` of the nonzero terms of the polynomial sorted in the decreasing order.

    monomials(vars::Tuple, degs::AbstractVector{Int}, filter::Function = m -> true)

Builds the vector of all the monovec `m` with variables `vars` such that the degree `deg(m)` is in `degs` and `filter(m)` is `true`.

### Examples

Calling `monomials` on ``4x^2y + xy + 2x`` should return an iterator of ``[x^2y, xy, x]``.

Calling `monomials((x, y), [1, 3], m -> exponent(m, y) != 1)` should return `[x^3, x*y^2, y^3, x]` where `x^2*y` and `y` have been excluded by the filter.
"""
monomials(p::APL) = monomial.(terms(p))

#$(SIGNATURES)
"""
    mindeg(p::AbstractPolynomialLike)

Returns the minimal total degree of the monomials of `p`, i.e. `minimum(deg, terms(p))`.

### Examples
Calling `mindeg` on on ``4x^2y + xy + 2x`` should return 1.
"""
function mindeg(p::AbstractPolynomialLike)
    minimum(deg, terms(p))
end

#$(SIGNATURES)
"""
    maxdeg(p::AbstractPolynomialLike)

Returns the maximal total degree of the monomials of `p`, i.e. `maximum(deg, terms(p))`.

### Examples
Calling `maxdeg` on on ``4x^2y + xy + 2x`` should return 3.
"""
function maxdeg(p::AbstractPolynomialLike)
    maximum(deg, terms(p))
end

#$(SIGNATURES)
"""
    extdeg(p::AbstractPolynomialLike)

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

#$(SIGNATURES)
"""
    leadingcoefficient(p::AbstractPolynomialLike)

Returns the coefficient of the leading term of `p`, i.e. `coefficient(leadingterm(p))`.

### Examples

Calling `leadingcoefficient` on ``4x^2y + xy + 2x`` should return ``4`` and calling it on ``0`` should return ``0``.
"""
function leadingcoefficient(p::AbstractPolynomialLike)
    coefficient(leadingterm(p))
end

#$(SIGNATURES)
"""
    leadingmonomial(p::AbstractPolynomialLike)

Returns the monomial of the leading term of `p`, i.e. `monomial(leadingterm(p))` or `first(monomials(p))`.

### Examples

Calling `leadingmonomial` on ``4x^2y + xy + 2x`` should return ``x^2y``.
"""
function leadingmonomial(p::AbstractPolynomialLike)
    # first(monomials(p)) would be more efficient for DynamicPolynomials but
    # monomial(leadingterm(p)) is more efficient for TypedPolynomials and is better if p is a term
    monomial(leadingterm(p))
end

#$(SIGNATURES)
"""
    removeleadingterm(p::AbstractPolynomialLike)

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
function removemonomials(p::AbstractPolynomialLike, mv::AbstractVector{MT}) where {MT <: AbstractMonomialLike}
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
    variables(p::AbstractPolynomialLike)

Returns the tuple of the variables of `p` in decreasing order. It could contain variables of zero degree, see the example section.

### Examples

Calling `variables(x^2*y)` should return `(x, y)` and calling `variables(x)` should return `(x,)`.
Note that the variables of `m` does not necessarily have nonzero exponent.
For instance, `variables([x^2*y, y*z][1])` is usually `(x, y, z)` since the two monomials have been promoted to a common type.
"""
function variables end

variables(x::AbstractVariable) = (x,)

"""
    nvariables(p::AbstractPolynomialLike)

Returns the number of variables in `p`, i.e. `length(variables(p))`. It could be more than the number of variables with nonzero exponent (see [the Examples section of `variables`](@ref MultivariatePolynomials.variables)).

### Examples

Calling `nvariables(x^2*y)` should return at least 2 and calling `nvariables(x)` should return at least 1.
"""
function nvariables(p::AbstractPolynomialLike)
    length(variables(p))
end

nvariables(::AbstractVariable) = 1

export constantterm, term, termtype, zeroterm, coefficient, monomial, powers, exponents, exponent, deg, isconstant, divides

function Base.hash(t::AbstractTerm, u::UInt)
    if iszero(t)
        hash(0, u)
    elseif coefficient(t) == 1
        hash(monomial(t), u)
    else
        hash(monomial(t), hash(coefficient(t), u))
    end
end

Base.zero(::Type{TT}) where {T, TT<:AbstractTermLike{T}} = zero(T) * constantmonomial(TT)
Base.one(::Type{TT}) where {T, TT<:AbstractTermLike{T}} = one(T) * constantmonomial(TT)
Base.zero(t::AbstractTermLike{T}) where {T} = zero(T) * constantmonomial(t)
Base.one(t::AbstractTermLike{T}) where {T} = one(T) * constantmonomial(t)

monomial(m::AbstractMonomial) = m

"""
    constantterm(α, p::AbstractPolynomialLike)

Creates a constant term with coefficient α and the same variables as p.

    constantterm{PT<:AbstractPolynomialType}(α, ::Type{PT}

Creates a constant term of the term type of a polynomial of type `PT`.
"""
constantterm(α, p) = α * constantmonomial(p)

"""
    term(p::AbstractPolynomialLike)

Converts the polynomial `p` to a term.
When applied on a polynomial, it throws an error if it has more than one term.
When applied to a term, it is the identity and does not copy it.
When applied to a monomial, it create a term of type `AbstractTerm{Int}`.
"""
function term(p::APL)
    if nterms(p) == 0
        zero(termtype(p))
    elseif nterms(p) > 1
        error("A polynomial is more than one term cannot be converted to a term.")
    else
        leadingterm(p)
    end
end
term(t::AbstractTerm) = t
term(m::AbstractMonomialLike) = 1 * m

termtype(p::APL) = Base.promote_op(first ∘ terms, p)

"""
    coefficient(t::AbstractTermLike)

Returns the coefficient of the term `t`.

### Examples

Calling `coefficient` on ``4x^2y`` should return ``4``.
"""
coefficient(m::AbstractMonomialLike) = 1

"""
    monomial(t::AbstractTermLike)

Returns the monomial of the term `t`.

### Examples

Calling `monomial` on ``4x^2y`` should return ``x^2y``.
"""
function monomial end

powers(t::AbstractTermLike) = tuplezip(vars(t), exponents(t))

"""
    exponent(t::AbstractTermLike, var::AbstractVariable)

Returns the exponent of the variable `var` in the monomial of the term `t`.

### Examples

Calling `exponent(x^2*y, x)` should return 2 and calling `exponent(x^2*y, y)` should return 1.
"""
exponent(t::AbstractTerm, var::AbstractVariable) = exponent(monomial(t), var)
exponent(v::AbstractVariable, var::AbstractVariable) = (v == var ? 1 : 0)

"""
    exponents(t::AbstractTermLike)

Returns the exponent of the variables in the monomial of the term `t`.

### Examples

Calling `exponents(x^2*y)` should return `(2, 1)`.
"""
exponents(t::AbstractTerm) = exponents(monomial(t))
exponents(v::AbstractVariable) = (1,)

"""
    deg(t::AbstractTermLike)

Returns the *total degree* of the monomial of the term `t`, i.e. `sum(exponents(t))`.

### Examples

Calling `deg(x^2*y)` should return 3 which is ``2 + 1``.
"""
deg(t::AbstractTermLike) = sum(exponents(t))

"""
    isconstant(m::AbstractMonomialLike)

Returns whether the monomial
"""
isconstant(t::AbstractTermLike) = deg(t) == 0
isconstant(v::AbstractVariable) = false

"""
    divides(t1::AbstractTermLike, t2::AbstractTermLike)

Returns whether `monomial(t1)` divides `monomial(t2)`.

### Examples

Calling `divides(2x^2y, 3xy)` should return false because `x^2y` does not divide `xy` since `x` has a degree 2 in `x^2y` which is greater than the degree of `x` on `xy`.
However, calling `divides(3xy, 2x^2y)` should return true.
"""
function divides end

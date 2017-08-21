export constantterm, term, termtype, zeroterm, coefficient, monomial, powers, exponents, degree, isconstant, divides

function Base.hash(t::AbstractTerm, u::UInt)
    if iszero(t)
        hash(0, u)
    elseif coefficient(t) == 1
        hash(monomial(t), u)
    else
        hash(monomial(t), hash(coefficient(t), u))
    end
end

# zero should return a polynomial since it is often used to keep the result of a summation of terms.
# For example, Base.vecdot(x::Vector{<:AbstractTerm}, y:Vector{Int}) starts with `s = zero(dot(first(x), first(y)))` and then adds terms.
# We want `s` to start as a polynomial for this operation to be efficient.
#Base.zero(::Type{TT}) where {T, TT<:AbstractTermLike{T}} = zero(T) * constantmonomial(TT)
#Base.zero(t::AbstractTermLike{T}) where {T} = zero(T) * constantmonomial(t)
zeroterm(::Type{TT}) where {T, TT<:APL{T}} = zero(T) * constantmonomial(TT)
zeroterm(t::APL{T}) where {T} = zero(T) * constantmonomial(t)

Base.zero(::Type{TT}) where {T, TT<:AbstractTermLike{T}} = zero(polynomialtype(TT))
Base.zero(t::AbstractTermLike{T}) where {T} = zero(polynomialtype(t))
Base.one(::Type{TT}) where {T, TT<:AbstractTermLike{T}} = one(T) * constantmonomial(TT)
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
function term(p::APL{T}) where T
    if nterms(p) == 0
        zeroterm(p)
    elseif nterms(p) > 1
        error("A polynomial is more than one term cannot be converted to a term.")
    else
        leadingterm(p)
    end
end
term(t::AbstractTerm) = t
term(m::AbstractMonomialLike) = 1 * m

"""
    termtype(p::AbstractPolynomialLike)

Returns the type of the terms of `p`.

    termtype(::Type{PT}) where PT<:AbstractPolynomialLike

Returns the type of the terms of a polynomial of type `PT`.

    termtype(p::AbstractPolynomialLike, ::Type{T}) where T

Returns the type of the terms of `p` but with coefficient type `T`.

    termtype(::Type{PT}, ::Type{T}) where {PT<:AbstractPolynomialLike, T}

Returns the type of the terms of a polynomial of type `PT` but with coefficient type `T`.
"""
termtype(p::P) where P <: APL = Base.promote_op(first ∘ terms, P)
termtype(::Type{V}, ::Type{C}) where {V <: AbstractVariable, C} = termtype(monomialtype(V), C)
termtype(::Type{T}) where T <: AbstractTerm = T
termtype(::Type{M}) where M<:AbstractMonomialLike = termtype(M, Int)

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

powers(t::AbstractTermLike) = tuplezip(variables(t), exponents(t))

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
_deg(v::AbstractVariable) = 0
_deg(v::AbstractVariable, power, powers...) = v == power[1] ? power[2] : _deg(v, powers...)
degree(t::AbstractTermLike, v::AbstractVariable) = _deg(v, powers(t)...)

degree(t::AbstractTerm, var::AbstractVariable) = degree(monomial(t), var)
degree(v::AbstractVariable, var::AbstractVariable) = (v == var ? 1 : 0)
function degree(m::AbstractMonomial, v::AbstractVariable)
    i = findfirst(variables(m), v)
    if iszero(i)
        0
    else
        exponents(m)[i]
    end
end


"""
    isconstant(m::AbstractMonomialLike)

Returns whether the monomial is constant.
"""
isconstant(t::AbstractTermLike) = all(iszero.(exponents(t)))
isconstant(v::AbstractVariable) = false

"""
    divides(t1::AbstractTermLike, t2::AbstractTermLike)

Returns whether `monomial(t1)` divides `monomial(t2)`.

### Examples

Calling `divides(2x^2y, 3xy)` should return false because `x^2y` does not divide `xy` since `x` has a degree 2 in `x^2y` which is greater than the degree of `x` on `xy`.
However, calling `divides(3xy, 2x^2y)` should return true.
"""
function divides end

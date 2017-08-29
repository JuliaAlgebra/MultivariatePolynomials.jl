export constantterm, term, termtype, zeroterm, coefficient, monomial

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
termtype(::Union{T, Type{T}}) where T <: AbstractTerm = T
termtype(::Union{P, Type{P}}) where P <: APL = Base.promote_op(first ∘ terms, P)
termtype(p::Union{APL, Type{<:APL}}, ::Type{T}) where T = termtype(termtype(p), T)
termtype(m::Union{M, Type{M}}) where M<:AbstractMonomialLike = termtype(m, Int)
termtype(v::Union{AbstractVariable, Type{<:AbstractVariable}}) = termtype(monomialtype(v))
termtype(v::Union{AbstractVariable, Type{<:AbstractVariable}}, ::Type{T}) where T = termtype(monomialtype(v), T)
termtype(::Union{AbstractVector{PT}, Type{<:AbstractVector{PT}}}) where PT <: APL = termtype(PT)
termtype(::Union{AbstractVector{PT}, Type{<:AbstractVector{PT}}}, ::Type{T}) where {PT <: APL, T} = termtype(PT, T)

"""
    coefficient(t::AbstractTermLike)

Returns the coefficient of the term `t`.

### Examples

Calling `coefficient` on ``4x^2y`` should return ``4``.
"""
coefficient(m::AbstractMonomialLike) = 1

"""
    coefficient(p::AbstractPolynomialLike)

Returns the coefficient type of `p`.

    coefficient(::Type{PT}) where PT

Returns the coefficient type of a polynomial of type `PT`.

### Examples

Calling `coefficienttype` on ``(4//5)x^2y`` should return `Rational{Int}`,
calling `coefficienttype` on ``1.0x^2y + 2.0x`` should return `Float64` and
calling `coefficienttype` on ``xy`` should return `Int`.
"""
coefficienttype(::Type{<:APL{T}}) where {T} = T
coefficienttype(::APL{T}) where {T} = T
#coefficienttype(::Type{T}) where {T} = T
#coefficienttype(::T) where {T} = T

"""
    monomial(t::AbstractTermLike)

Returns the monomial of the term `t`.

### Examples

Calling `monomial` on ``4x^2y`` should return ``x^2y``.
"""
function monomial end
monomial(m::AbstractMonomial) = m

"""
    constantterm(α, p::AbstractPolynomialLike)

Creates a constant term with coefficient α and the same variables as p.

    constantterm(α, ::Type{PT} where {PT<:AbstractPolynomialType}

Creates a constant term of the term type of a polynomial of type `PT`.
"""
constantterm(α, p) = α * constantmonomial(p)

# zero should return a polynomial since it is often used to keep the result of a summation of terms.
# For example, Base.vecdot(x::Vector{<:AbstractTerm}, y:Vector{Int}) starts with `s = zero(dot(first(x), first(y)))` and then adds terms.
# We want `s` to start as a polynomial for this operation to be efficient.
#Base.zero(::Type{TT}) where {T, TT<:AbstractTermLike{T}} = zero(T) * constantmonomial(TT)
#Base.zero(t::AbstractTermLike{T}) where {T} = zero(T) * constantmonomial(t)
"""
    zeroterm(p::AbstractPolynomialLike{T}) where T

Equivalent to `constantterm(zero(T), p)`.

    zeroterm(α, ::Type{PT} where {T, PT<:AbstractPolynomialType{T}}

Equivalent to `constantterm(zero(T), PT)`.
"""
zeroterm(::Type{PT}) where {T, PT<:APL{T}} = constantterm(zero(T), PT)
zeroterm(p::APL{T}) where {T} = constantterm(zero(T), p)

Base.zero(::Type{TT}) where {T, TT<:AbstractTermLike{T}} = zero(polynomialtype(TT))
Base.zero(t::AbstractTermLike{T}) where {T} = zero(polynomialtype(t))
Base.one(::Type{TT}) where {T, TT<:AbstractTermLike{T}} = one(T) * constantmonomial(TT)
Base.one(t::AbstractTermLike{T}) where {T} = one(T) * constantmonomial(t)

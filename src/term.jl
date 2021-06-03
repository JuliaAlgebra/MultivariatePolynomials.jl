export constantterm, term, termtype, zeroterm, coefficient, monomial

struct Term{T, M<:AbstractMonomial} <: AbstractTerm{T}
    coefficient::T
    monomial::M
end

coefficient(t::Term) = t.coefficient
monomial(t::Term) = t.monomial
termtype(::Union{Term{C, M}, Type{Term{C, M}}}, ::Type{T}) where {C, M, T} = Term{T, M}
termtype(::Union{M, Type{M}}, ::Type{T}) where {M<:AbstractMonomial, T} = Term{T, M}
monomialtype(::Union{Term{T, M}, Type{Term{T, M}}}) where {T, M} = M

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
When applied to a monomial, it create a term of type `AbstractTerm{Int}`.
"""
function term end
term(coef, mono::AbstractMonomial) = Term(coef, mono)
term(coef, var::AbstractVariable) = term(coef, monomial(var))
term(p::APL) = convert(termtype(p), p)

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
function coefficient(p::AbstractPolynomialLike{T}, m::AbstractMonomialLike) where T
    for t in terms(p)
        if monomial(t) == m
            return coefficient(t)
        end
    end
    zero(T)
end

"""
    coefficient(p::AbstractPolynomialLike, m::AbstractMonomialLike, vars)::AbstractPolynomialLike

Returns the coefficient of the monomial `m` of the polynomial `p` considered as a polynomial in variables
`vars`.

### Example
Calling `coefficient((a+b)x^2+2x+y*x^2, x^2, [x,y])` should return `a+b`.
Calling `coefficient((a+b)x^2+2x+y*x^2, x^2, [x])` should return `a+b+y`.
"""
function coefficient(f::APL, m::AbstractMonomialLike, vars)
    coeff = zero(f)
    for t in terms(f)
        match = true
        for v in vars
            if degree(t, v) != degree(m, v)
                match = false
                break
            end
        end
        match || continue

        coeff += coefficient(t) * (_div(monomial(t), m))
    end
    coeff
end

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
coefficienttype(::Union{PT, Type{PT}, AbstractVector{PT}, Type{<:AbstractVector{PT}}}) where {T, PT<:APL{T}} = T
#coefficienttype(::{T, Type{T}}) where {T} = T

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
    constantterm(α, p::AbstractPolynomialLike)

Creates a constant term with coefficient α and the same variables as p.

    constantterm(α, ::Type{PT} where {PT<:AbstractPolynomialLike}

Creates a constant term of the term type of a polynomial of type `PT`.
"""
constantterm(α, p) = term(α, constantmonomial(p))

# zero should return a polynomial since it is often used to keep the result of a summation of terms.
# For example, Base.vecdot(x::Vector{<:AbstractTerm}, y:Vector{Int}) starts with `s = zero(dot(first(x), first(y)))` and then adds terms.
# We want `s` to start as a polynomial for this operation to be efficient.
#Base.zero(::Type{TT}) where {T, TT<:AbstractTermLike{T}} = zero(T) * constantmonomial(TT)
#Base.zero(t::AbstractTermLike{T}) where {T} = zero(T) * constantmonomial(t)
"""
    zeroterm(p::AbstractPolynomialLike{T}) where T

Equivalent to `constantterm(zero(T), p)`.

    zeroterm(α, ::Type{PT} where {T, PT<:AbstractPolynomialLike{T}}

Equivalent to `constantterm(zero(T), PT)`.
"""
zeroterm(::Type{PT}) where {T, PT<:APL{T}} = constantterm(zero(T), PT)
zeroterm(p::APL{T}) where {T} = constantterm(zero(T), p)

Base.zero(::Type{TT}) where {T, TT<:AbstractTermLike{T}} = zero(polynomialtype(TT))
Base.zero(t::AbstractTermLike{T}) where {T} = zero(polynomialtype(t))
MA.promote_operation(::typeof(zero), PT::Type{<:AbstractTermLike}) = polynomialtype(PT)
Base.one(::Type{TT}) where {T, TT<:AbstractTermLike{T}} = one(T) * constantmonomial(TT)
Base.one(t::AbstractTermLike{T}) where {T} = one(T) * constantmonomial(t)
MA.promote_operation(::typeof(one), TT::Type{<:AbstractTermLike}) = termtype(TT)
MA.promote_operation(::typeof(one), PT::Type{<:AbstractPolynomialLike}) = polynomialtype(PT)

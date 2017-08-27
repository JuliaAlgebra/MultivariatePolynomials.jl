export polynomial, polynomialtype, terms, nterms, coefficients, monomials
export coefficienttype, monomialtype
export mindegree, maxdegree, extdegree
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

coefficienttype(::Type{<:APL{T}}) where {T} = T
coefficienttype(::APL{T}) where {T} = T
#coefficienttype(::Type{T}) where {T} = T
#coefficienttype(::T) where {T} = T

changecoefficienttype(::Type{TT}, ::Type{T}) where {TT<:AbstractTermLike, T} = termtype(TT, T)
changecoefficienttype(::Type{PT}, ::Type{T}) where {PT<:AbstractPolynomial, T} = polynomialtype(PT, T)

changecoefficienttype(p::PT, ::Type{T}) where {PT<:APL, T} = convert(changecoefficienttype(PT, T), p)

abstract type ListState end
abstract type UnsortedState <: ListState end
struct MessyState <: UnsortedState end
# No duplicates or zeros
struct UniqState <: UnsortedState end
sortstate(::MessyState) = SortedState()
sortstate(::UniqState) = SortedUniqState()
struct SortedState <: ListState end
struct SortedUniqState <: ListState end

"""
    polynomial(p::AbstractPolynomialLike)

Converts `p` to a value with polynomial type.

    polynomial(p::AbstractPolynomialLike, ::Type{T}) where T

Converts `p` to a value with polynomial type with coefficient type `T`.

    polynomial(a::AbstractVector, mv::AbstractVector{<:AbstractMonomialLike})

Creates a polynomial equal to `dot(a, mv)`.

    polynomial(terms::AbstractVector{<:AbstractTerm}, s::ListState=MessyState())

Creates a polynomial equal to `sum(terms)` where `terms` are guaranteed to be in state `s`.

    polynomial(f::Function, mv::AbstractVector{<:AbstractMonomialLike})

Creates a polynomial equal to `sum(f(i) * mv[i] for i in 1:length(mv))`.

### Examples

Calling `polynomial([2, 4, 1], [x, x^2*y, x*y])` should return ``4x^2y + xy + 2x``.
"""
polynomial(p::AbstractPolynomial) = p
polynomial(p::APL{T}, ::Type{T}) where T = polynomial(terms(p))
polynomial(p::APL{T}) where T = polynomial(p, T)
polynomial(ts::AbstractVector, s::ListState=MessyState()) = sum(ts)
polynomial(ts::AbstractVector{<:AbstractTerm}, s::SortedUniqState) = polynomial(coefficient.(ts), monomial.(ts), s)
polynomial(a::AbstractVector, x::AbstractVector, s::ListState=MessyState()) = polynomial([α * m for (α, m) in zip(a, x)], s)
polynomial(f::Function, mv::AbstractVector{<:AbstractMonomialLike}) = polynomial([f(i) * mv[i] for i in 1:length(mv)])
function polynomial(Q::AbstractMatrix, mv::AbstractVector)
    dot(mv, Q * mv)
end
function polynomial(Q::AbstractMatrix, mv::AbstractVector, ::Type{T}) where T
    polynomial(polynomial(Q, mv), T)
end

polynomialtype(::Union{P, Type{P}}) where P<:APL = Base.promote_op(polynomial, P)
polynomialtype(::Type{P}) where P<:AbstractPolynomial = P
polynomialtype(::Type{M}, ::Type{T}) where {M<:AbstractMonomialLike, T} = polynomialtype(termtype(M, T))
polynomialtype(::Type{P}, ::Type{T}) where {P, T} = polynomialtype(polynomialtype(P), T)
polynomialtype(p, ::Type{T}) where T = polynomialtype(typeof(p), T)

function uniqterms(ts::AbstractVector{T}) where T <: AbstractTerm
    result = T[]
    sizehint!(result, length(ts))
    for t in ts
        if !iszero(t)
            if isempty(result) || monomial(t) != monomial(last(result))
                push!(result, t)
            else
                coef = coefficient(last(result)) + coefficient(t)
                if iszero(coef)
                    pop!(result)
                else
                    result[end] = coef * monomial(t)
                end
            end
        end
    end
    result
end
polynomial(ts::AbstractVector{<:AbstractTerm}, s::SortedState) = polynomial(uniqterms(ts), SortedUniqState())
polynomial(ts::AbstractVector{<:AbstractTerm}, s::UnsortedState=MessyState()) = polynomial(sort(ts, lt=(>)), sortstate(s))

"""
    terms(p::AbstractPolynomialLike)

Returns an iterator over the nonzero terms of the polynomial `p` sorted in the decreasing monomial order.

### Examples

Calling `terms` on ``4x^2y + xy + 2x`` should return an iterator of ``[4x^2y, xy, 2x]``.
"""
terms(t::AbstractTermLike) = iszero(t) ? termtype(t)[] : [term(t)]
terms(p::AbstractPolynomialLike) = terms(polynomial(p))

"""
    nterms(p::AbstractPolynomialLike)

Returns the number of nonzero terms in `p`, i.e. `length(terms(p))`.

### Examples

Calling `nterms` on ``4x^2y + xy + 2x`` should return 3.
"""
function nterms end

nterms(t::AbstractTermLike) = iszero(t) ? 0 : 1
nterms(p::AbstractPolynomialLike) = length(terms(p))

"""
    coefficients(p::AbstractPolynomialLike)

Returns an iterator over the coefficients of `p` of the nonzero terms of the polynomial sorted in the decreasing monomial order.

    coefficients(p::AbstractPolynomialLike, X::AbstractVector)

Returns an iterator over the coefficients of the monomials of `X` in `p` where `X` is a monomial vector not necessarily sorted but with no duplicate entry.

### Examples

Calling `coefficients` on ``4x^2y + xy + 2x`` should return an iterator of ``[4, 1, 2]``.
Calling `coefficients(4x^2*y + x*y + 2x + 3, [x, 1, x*y, y])` should return an iterator of ``[2, 3, 1, 0]``.
"""
coefficients(p::APL) = coefficient.(terms(p))
function coefficients(p::APL{T}, X::AbstractVector) where T
    σ, mv = sortmonovec(X)
    @assert length(mv) == length(X) # no duplicate in X
    c = Vector{T}(length(mv))
    i = 1
    for t in terms(p)
        m = monomial(t)
        while i <= length(mv) && mv[i] > m
            c[σ[i]] = zero(T)
            i += 1
        end
        if i <= length(mv) && mv[i] == m
            c[σ[i]] = coefficient(t)
            i += 1
        end
    end
    c
end

"""
    monomials(p::AbstractPolynomialLike)

Returns an iterator over the monomials of `p` of the nonzero terms of the polynomial sorted in the decreasing order.

    monomials(vars::Tuple, degs::AbstractVector{Int}, filter::Function = m -> true)

Builds the vector of all the monovec `m` with variables `vars` such that the degree `degree(m)` is in `degs` and `filter(m)` is `true`.

### Examples

Calling `monomials` on ``4x^2y + xy + 2x`` should return an iterator of ``[x^2y, xy, x]``.

Calling `monomials((x, y), [1, 3], m -> degree(m, y) != 1)` should return `[x^3, x*y^2, y^3, x]` where `x^2*y` and `y` have been excluded by the filter.
"""
monomials(p::APL) = monomial.(terms(p))

#$(SIGNATURES)
"""
    mindegree(p::Union{AbstractPolynomialLike, AbstractVector{<:AbstractTermLike}})

Returns the minimal total degree of the monomials of `p`, i.e. `minimum(degree, terms(p))`.

    mindegree(p::Union{AbstractPolynomialLike, AbstractVector{<:AbstractTermLike}}, v::AbstractVariable)

Returns the minimal degree of the monomials of `p` in the variable `v`, i.e. `minimum(degree.(terms(p), v))`.

### Examples
Calling `mindegree` on on ``4x^2y + xy + 2x`` should return 1, `mindegree(4x^2y + xy + 2x, x)` should return 1 and  `mindegree(4x^2y + xy + 2x, y)` should return 0.
"""
function mindegree(X::AbstractVector{<:AbstractTermLike}, args...)
    minimum(t -> degree(t, args...), X)
end
function mindegree(p::AbstractPolynomialLike, args...)
    mindegree(terms(p), args...)
end

#$(SIGNATURES)
"""
    maxdegree(p::Union{AbstractPolynomialLike, AbstractVector{<:AbstractTermLike}})

Returns the maximal total degree of the monomials of `p`, i.e. `maximum(degree, terms(p))`.

    maxdegree(p::Union{AbstractPolynomialLike, AbstractVector{<:AbstractTermLike}}, v::AbstractVariable)

Returns the maximal degree of the monomials of `p` in the variable `v`, i.e. `maximum(degree.(terms(p), v))`.

### Examples
Calling `maxdegree` on on ``4x^2y + xy + 2x`` should return 3, `maxdegree(4x^2y + xy + 2x, x)` should return 2 and  `maxdegree(4x^2y + xy + 2x, y)` should return 1.
"""
function maxdegree(X::AbstractVector{<:AbstractTermLike}, args...)
    maximum(t -> degree(t, args...), X)
end
function maxdegree(p::AbstractPolynomialLike, args...)
    maxdegree(terms(p), args...)
end

#$(SIGNATURES)
"""
    extdegree(p::Union{AbstractPolynomialLike, AbstractVector{<:AbstractTermLike}})

Returns the extremal total degrees of the monomials of `p`, i.e. `(mindegree(p), maxdegree(p))`.

    extdegree(p::Union{AbstractPolynomialLike, AbstractVector{<:AbstractTermLike}}, v::AbstractVariable)

Returns the extremal degrees of the monomials of `p` in the variable `v`, i.e. `(mindegree(p, v), maxdegree(p, v))`.

### Examples
Calling `extdegree` on on ``4x^2y + xy + 2x`` should return `(1, 3)`, `extdegree(4x^2y + xy + 2x, x)` should return `(1, 2)` and  `maxdegree(4x^2y + xy + 2x, y)` should return `(0, 1)`.
"""
function extdegree(p::Union{AbstractPolynomialLike, AbstractVector{<:AbstractTermLike}}, args...)
    (mindegree(p, args...), maxdegree(p, args...))
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
    # Iterators.drop returns an Interators.Drop which is not an AbstractVector
    polynomial(terms(p)[2:end], SortedUniqState())
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
        if i > length(smv) || smv[i] != m
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
Note that the variables of `m` does not necessarily have nonzero degree.
For instance, `variables([x^2*y, y*z][1])` is usually `(x, y, z)` since the two monomials have been promoted to a common type.
"""
function variables end

"""
    nvariables(p::AbstractPolynomialLike)

Returns the number of variables in `p`, i.e. `length(variables(p))`. It could be more than the number of variables with nonzero degree (see [the Examples section of `variables`](@ref MultivariatePolynomials.variables)).

### Examples

Calling `nvariables(x^2*y)` should return at least 2 and calling `nvariables(x)` should return at least 1.
"""
nvariables(::Union{AbstractVariable, Type{<:AbstractVariable}}) = 1

export polynomial, polynomial!, polynomialtype, terms, nterms, coefficients, monomials
export coefficienttype, monomialtype
export mindegree, maxdegree, extdegree, effective_variables
export leadingterm, leadingcoefficient, leadingmonomial
export removeleadingterm, removemonomials, monic
export mapcoefficients, mapcoefficients!, mapcoefficients_to!

LinearAlgebra.norm(p::AbstractPolynomialLike, r::Int=2) = LinearAlgebra.norm(coefficients(p), r)

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
function polynomial end
function polynomial(p::APL, args::Vararg{Type,N}) where {N}
    return polynomial!(copy(p), args...)
end
function polynomial(Q::AbstractMatrix, mv::AbstractVector)
    LinearAlgebra.dot(mv, Q * mv)
end
function polynomial(Q::AbstractMatrix, mv::AbstractVector, ::Type{T}) where T
    polynomial(polynomial(Q, mv), T)
end

polynomial(f::F, mv::AbstractVector{<:AbstractMonomialLike}) where F = polynomial!([term(f(i), mv[i]) for i in 1:length(mv)])

function polynomial(a::AbstractVector, x::AbstractVector, s::ListState=MessyState())
    # If `x` is e.g. `[v, 1]` then it will contains terms that are convertible to monomials.
    return polynomial([term(α, convert(monomialtype(m), m)) for (α, m) in zip(a, x)], s)
end

polynomial(ts::AbstractVector, s::ListState=MessyState()) = sum(ts)

function polynomial!(p::APL, args::Vararg{Any,N}) where N
    return convert(polynomialtype(p, args...), p)
end

polynomial!(ts::AbstractVector, s::ListState=MessyState()) = sum(ts)

function polynomial!(ts::AbstractVector{TT}, s::SortedUniqState) where {TT<:AbstractTerm}
    return polynomialtype(TT)(ts)
end

function polynomial!(ts::AbstractVector{<:AbstractTerm}, s::SortedState)
    polynomial!(uniqterms!(ts), SortedUniqState())
end
function polynomial!(ts::AbstractVector{<:AbstractTerm}, s::UnsortedState=MessyState())
    polynomial!(sort!(ts, lt=(>)), sortstate(s))
end

_collect(v::Vector) = v
_collect(v::AbstractVector) = collect(v)
polynomial(ts::AbstractVector{<:AbstractTerm}, args::Vararg{ListState, N}) where {N} = polynomial!(MA.mutable_copy(_collect(ts)), args...)

"""
    polynomialtype(p::AbstractPolynomialLike)

Returns the type that `p` would have if it was converted into a polynomial.

    polynomialtype(::Type{PT}) where PT<:AbstractPolynomialLike

Returns the same as `polynomialtype(::PT)`.

    polynomialtype(p::AbstractPolynomialLike, ::Type{T}) where T

Returns the type that `p` would have if it was converted into a polynomial of coefficient type `T`.

    polynomialtype(::Type{PT}, ::Type{T}) where {PT<:AbstractPolynomialLike, T}

Returns the same as `polynomialtype(::PT, ::Type{T})`.
"""
function polynomialtype end
polynomialtype(::Type{T}) where T <: AbstractTerm = error("`polynomialtype` not implemented for $T")
polynomialtype(::Union{P, Type{P}}) where P <: APL = polynomialtype(termtype(P))
polynomialtype(::Union{P, Type{P}}) where P <: AbstractPolynomial = P
polynomialtype(::Union{M, Type{M}}) where M<:AbstractMonomialLike = polynomialtype(termtype(M))
polynomialtype(::Union{M, Type{M}}, ::Type{T}) where {M<:AbstractMonomialLike, T} = polynomialtype(termtype(M, T))
polynomialtype(::Union{P, Type{P}}, ::Type{T}) where {P <: APL, T} = polynomialtype(polynomialtype(P), T)
polynomialtype(::Union{AbstractVector{PT}, Type{<:AbstractVector{PT}}}) where PT <: APL = polynomialtype(PT)
polynomialtype(::Union{AbstractVector{PT}, Type{<:AbstractVector{PT}}}, ::Type{T}) where {PT <: APL, T} = polynomialtype(PT, T)

function uniqterms!(ts::AbstractVector{<:AbstractTerm})
    i = firstindex(ts)
    for j in Iterators.drop(eachindex(ts), 1)
        if !iszero(ts[j])
            if monomial(ts[i]) == monomial(ts[j])
                ts[i] = term(MA.add!!(coefficient(ts[i]), coefficient(ts[j])), monomial(ts[i]))
            else
                if !iszero(ts[i])
                    i += 1
                end
                ts[i] = MA.copy_if_mutable(ts[j])
            end
        end
    end
    if i < length(ts)
        if iszero(ts[i])
            i -= 1
        end
        resize!(ts, i)
    end
    ts
end

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
    c = zeros(T, length(mv))
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
monomials(p::APL) = monovec(monomial.(terms(p)))

function isconstant(p::APL)
    n = nterms(p)
    return iszero(n) || (isone(n) && isconstant(first(terms(p))))
end

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
    isempty(X) ? 0 : minimum(t -> degree(t, args...), X)
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
function maxdegree(X::AbstractVector{<:AbstractTermLike}, args::Vararg{Any,N}) where {N}
    return mapreduce(t -> degree(t, args...), max, X, init=0)
end
function maxdegree(p::AbstractPolynomialLike, args::Vararg{Any,N}) where {N}
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
    effective_variables(p::AbstractPolynomialLike)

Return a vector of `eltype` `variable_union_type(p)` (see [`variable_union_type`](@ref)),
containing all the variables that has nonzero degree in at least one term.
That is, return all the variables `v` such that `maxdegree(p, v)` is not zero.
The returned vector is sorted in decreasing order.
"""
function effective_variables(p::AbstractPolynomialLike)
    VT = variable_union_type(p)
    return VT[v for v in variables(p) if !iszero(maxdegree(p, v))]
end

"""
    leadingterm(p::AbstractPolynomialLike)

Returns the coefficient of the leading term, i.e. `first(terms(p))`.

### Examples

Calling `leadingterm` on ``4x^2y + xy + 2x`` should return ``4x^2y``.
"""
function leadingterm(p::AbstractPolynomialLike)
    if iszero(p)
        zeroterm(p)
    else
        first(terms(p))
    end
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
function MA.promote_operation(::typeof(removeleadingterm), ::Type{PT}) where {PT<:AbstractPolynomial}
    return PT
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
    monic(p::AbstractPolynomialLike)

Returns `p / leadingcoefficient(p)` where the leading coefficient of the returned polynomials is made sure to be exactly one to avoid rounding error.
"""
function monic(p::APL)
    α = leadingcoefficient(p)
    polynomial!(_divtoone.(terms(p), α))
end
monic(m::AbstractMonomialLike) = m
monic(t::AbstractTermLike{T}) where T = term(one(T), monomial(t))

function _divtoone(t::AbstractTermLike{T}, α::S) where {T, S}
    U = Base.promote_op(/, T, S)
    β = coefficient(t)
    if β == α
        term(one(U), monomial(t))
    else
        term((β / α), monomial(t))
    end
end

# TODO deprecate
mapcoefficientsnz(f::F, p::APL) where F = mapcoefficients(f, p, nonzero = true)
mapcoefficientsnz_to!(output::APL, f::F, p::APL) where F = mapcoefficients_to!(output, f, p, nonzero = true)

"""
    mapcoefficients(f::Function, p::AbstractPolynomialLike, nonzero = false)

Returns a polynomial with the same monomials as `p` but each coefficient `α` is replaced by `f(α)`.
The function may return zero in which case the term is dropped.
If the function is known to never returns zero for a nonzero input, `nonzero`
can be set to `true` to get a small speedup.

See also [`mapcoefficients!`](@ref) and [`mapcoefficients_to!`](@ref).

### Examples

Calling `mapcoefficients(α -> mod(3α, 6), 2x*y + 3x + 1)` should return `3x + 3`.
"""
function mapcoefficients end
function mapcoefficients(f::F, p::AbstractPolynomialLike; nonzero = false) where F # Not used by either TypedPolynomials or DynamicPolynomials but used by CustomPoly in tests. FIXME Remove in a breaking release
    # Invariant: p has only nonzero coefficient
    # therefore f(α) will be nonzero for every coefficient α of p
    # hence we can use Uniq
    polynomial!(mapcoefficients.(f, terms(p)), nonzero ? SortedUniqState() : SortedState())
end
function mapcoefficients(f::F, t::AbstractTermLike; nonzero = false) where F
    return term(f(coefficient(t)), monomial(t))
end

"""
    mapcoefficients!(f::Function, p::AbstractPolynomialLike, nonzero = false)

Mutate `p` by replacing each coefficient `α` by `f(α)`.
The function may return zero in which case the term is dropped.
If the function is known to never returns zero for a nonzero input, `nonzero`
can be set to `true` to get a small speedup.
The function returns `p`, which is identically equal to the second argument.

See also [`mapcoefficients`](@ref) and [`mapcoefficients_to!`](@ref).

### Examples

Let `p = 2x*y + 3x + 1`, after `mapcoefficients!(α -> mod(3α, 6), p)`, `p` is
equal to `3x + 3`.
"""
function mapcoefficients! end

"""
    mapcoefficients_to!(output::AbstractPolynomialLike, f::Function, p::AbstractPolynomialLike, nonzero = false)

Mutate `output` by replacing each coefficient `α` of `p` by `f(α)`.
The function may return zero in which case the term is dropped.
If the function is known to never returns zero for a nonzero input, `nonzero`
can be set to `true` to get a small speedup.
The function returns `output`, which is identically equal to the first argument.

See also [`mapcoefficients!`](@ref) and [`mapcoefficients`](@ref).
"""
function mapcoefficients_to! end

"""
    deg_num_leading_terms(p::AbstractPolynomialLike, var)

Return `deg, num` where `deg = maxdegree(p, var)` and `num` is the number of
terms `t` such that `degree(t, var) == deg`.
"""
function deg_num_leading_terms(p::AbstractPolynomialLike, var)
    deg = 0
    num = 0
    for mono in monomials(p)
        d = degree(mono, var)
        if d > deg
            deg = d
            num = 1
        elseif d == deg
            num += 1
        end
    end
    return deg, num
end

Base.ndims(::Union{Type{<:AbstractPolynomialLike}, AbstractPolynomialLike}) = 0
Base.broadcastable(p::AbstractPolynomialLike) = Ref(p)

Base.iszero(v::AbstractVariable) = false
Base.iszero(m::AbstractMonomial) = false
Base.iszero(t::AbstractTerm) = iszero(coefficient(t))
Base.iszero(t::AbstractPolynomial) = iszero(nterms(t))

Base.isone(v::AbstractVariable) = false
Base.isone(m::AbstractMonomial) = isconstant(m)
Base.isone(t::AbstractTerm) = isone(coefficient(t)) && isconstant(monomial(t))
function Base.isone(p::AbstractPolynomial)
    return isone(nterms(p)) && isone(first(terms(p)))
end

Base.:(==)(α::Dict, x::_APL) = false
Base.:(==)(x::_APL, α::Dict) = false
Base.:(==)(α::Dict, x::RationalPoly) = false
Base.:(==)(x::RationalPoly, α::Dict) = false

function right_term_eq(p::AbstractPolynomial, t; comp = (==))
    if iszero(p)
        iszero(t)
    else
        # terms/nterms ignore zero terms
        nterms(p) == 1 && comp(leading_term(p), t)
    end
end
right_term_eq(p::_APL, t; comp = (==)) = right_term_eq(polynomial(p), t; comp)

left_constant_eq(α, v::AbstractVariable; comp = (==)) = false
right_constant_eq(v::AbstractVariable, α; comp = (==)) = false
function _term_constant_eq(t::AbstractTermLike, α; comp = (==))
    if iszero(t)
        iszero(α)
    else
        comp(α, coefficient(t)) && isconstant(t)
    end
end
function left_constant_eq(α, t::AbstractTermLike; comp = (==))
    return _term_constant_eq(t, α; comp)
end
function right_constant_eq(t::AbstractTermLike, α; comp = (==))
    return _term_constant_eq(t, α; comp)
end
left_constant_eq(α, p::_APL; comp = (==)) = right_term_eq(p, α; comp)
right_constant_eq(p::_APL, α; comp = (==)) = right_term_eq(p, α; comp)

function Base.:(==)(mono::AbstractMonomial, v::AbstractVariable)
    return isone(degree(mono)) && variable(mono) == v
end
function Base.:(==)(v::AbstractVariable, mono::AbstractMonomial)
    return isone(degree(mono)) && v == variable(mono)
end
function Base.:(==)(t::AbstractTerm, mono::AbstractMonomialLike)
    return isone(coefficient(t)) && monomial(t) == mono
end
function Base.:(==)(mono::AbstractMonomialLike, t::AbstractTerm)
    return isone(coefficient(t)) && mono == monomial(t)
end

function _compare_term(t1::AbstractTerm, t2::AbstractTerm, comp)
    c1 = coefficient(t1)
    c2 = coefficient(t2)
    if iszero(c1)
        iszero(c2)
    else
        comp(c1, c2) && comp(monomial(t1), monomial(t2))
    end
end

Base.:(==)(t1::AbstractTerm, t2::AbstractTerm) = _compare_term(t1, t2, ==)
Base.:(==)(p::AbstractPolynomial, t::AbstractTerm) = right_term_eq(p, t)
Base.:(==)(t::AbstractTerm, p::AbstractPolynomial) = right_term_eq(p, t)
function Base.isequal(t1::AbstractTerm, t2::AbstractTerm)
    return _compare_term(t1, t2, isequal)
end
function Base.isequal(p::AbstractPolynomial, t::AbstractTerm)
    return right_term_eq(p, t; comp = isequal)
end
function Base.isequal(t::AbstractTerm, p::AbstractPolynomial)
    return right_term_eq(p, t; comp = isequal)
end

function compare_terms(p1::AbstractPolynomial, p2::AbstractPolynomial, isz, op)
    i1 = 1
    i2 = 1
    t1 = terms(p1)
    t2 = terms(p2)
    while true
        while i1 <= length(t1) && isz(coefficient(t1[i1]))
            i1 += 1
        end
        while i2 <= length(t2) && isz(coefficient(t2[i2]))
            i2 += 1
        end
        if i1 > length(t1) && i2 > length(t2)
            return true
        end
        if i1 > length(t1) || i2 > length(t2)
            return false
        end
        if !op(t1[i1], t2[i2])
            return false
        end
        i1 += 1
        i2 += 1
    end
end

# Can there be zero term in TypedPolynomials ?
#function (==)(p1::AbstractPolynomial, p2::AbstractPolynomial)
#    nterms(p1) != nterms(p2) && return false
#    for (t1, t2) in zip(terms(p1), terms(p2))
#        @assert !iszero(t1) && !iszero(t2) # There should be no zero term
#        if t1 != t2
#            return false
#        end
#    end
#    return true
#end
function Base.:(==)(p1::AbstractPolynomial, p2::AbstractPolynomial)
    return compare_terms(p1, p2, iszero, ==)
end
function Base.isequal(p1::AbstractPolynomial, p2::AbstractPolynomial)
    return compare_terms(p1, p2, iszero, isequal)
end

Base.:(==)(p::RationalPoly, q::RationalPoly) = p.num * q.den == q.num * p.den
# Solve ambiguity with (::PolyType, ::Any)
Base.:(==)(p::_APL, q::RationalPoly) = p * q.den == q.num
Base.:(==)(q::RationalPoly, p::_APL) = p == q
Base.:(==)(α::_Constant, q::RationalPoly) = α * q.den == q.num
Base.:(==)(q::RationalPoly, α::_Constant) = α == q
function Base.isequal(p::RationalPoly, q::RationalPoly)
    return isequal(p.num * q.den, q.num * p.den)
end
# Solve ambiguity with (::PolyType, ::Any)
Base.isequal(p::_APL, q::RationalPoly) = isequal(p * q.den, q.num)
Base.isequal(q::RationalPoly, p::_APL) = isequal(p, q)
Base.isequal(α::_Constant, q::RationalPoly) = isequal(α * q.den, q.num)
Base.isequal(q::RationalPoly, α::_Constant) = isequal(α, q)

# α could be a JuMP affine expression
isapproxzero(α; ztol::Real = 0.0) = false
function isapproxzero(α::Number; ztol::Real = Base.rtoldefault(α, α, 0))
    return abs(α) <= ztol
end

isapproxzero(m::AbstractMonomialLike; kwargs...) = false
function isapproxzero(t::AbstractTermLike; kwargs...)
    return isapproxzero(coefficient(t); kwargs...)
end
function isapproxzero(p::_APL; kwargs...)
    return all(term -> isapproxzero(term; kwargs...), terms(p))
end
isapproxzero(p::RationalPoly; kwargs...) = isapproxzero(p.num; kwargs...)

function Base.isapprox(t1::AbstractTerm, t2::AbstractTerm; kwargs...)
    return isapprox(coefficient(t1), coefficient(t2); kwargs...) &&
           monomial(t1) == monomial(t2)
end
function Base.isapprox(
    p1::AbstractPolynomial{S},
    p2::AbstractPolynomial{T};
    atol = 0,
    ztol::Real = iszero(atol) ? Base.rtoldefault(S, T, 0) : atol,
    kwargs...,
) where {S,T}
    return compare_terms(
        p1,
        p2,
        t -> isapproxzero(t; ztol = ztol),
        (x, y) -> isapprox(x, y; atol = atol, kwargs...),
    )
end

function Base.isapprox(p::RationalPoly, q::RationalPoly; kwargs...)
    return isapprox(p.num * q.den, q.num * p.den; kwargs...)
end
function Base.isapprox(p::RationalPoly, q::_APL; kwargs...)
    return isapprox(p.num, q * p.den; kwargs...)
end
function Base.isapprox(p::_APL, q::RationalPoly; kwargs...)
    return isapprox(p * q.den, q.num; kwargs...)
end
function Base.isapprox(q::RationalPoly{C}, α; kwargs...) where {C}
    return isapprox(q, constant_term(α, q.den); kwargs...)
end
function Base.isapprox(α, q::RationalPoly{C}; kwargs...) where {C}
    return isapprox(constant_term(α, q.den), q; kwargs...)
end

# TODO refer to the parallel with `mstructure(p)(e1, e2)` which gives the result
#      of multiplying the monomials corresponding to the exponent vectors `e1`
#      and `e2`.
"""
    abstract type AbstractMonomialOrdering end

Abstract type for monomial ordering as defined in [CLO13, Definition 2.2.1, p. 55]

Given an ordering `ordering::AbstractMonomialOrdering` and vector of exponents `e1`
and `e2`, `cmp(ordering, e1, e2)` returns a negative number if `e1` is before `e2`
in the ordering, a positive number if `e2` is before `e1` and 0 if they are equal.
For convenience, `ordering(e1, e2)` returns a `Bool` indicating whether
`cmp(ordering, e1, e2)` is negative.

[CLO13] Cox, D., Little, J., & OShea, D.
*Ideals, varieties, and algorithms: an introduction to computational algebraic geometry and commutative algebra*.
Springer Science & Business Media, **2013**.
"""
abstract type AbstractMonomialOrdering end

# We can't write this with a type instead of an instance so this motivates
# why we work with instances and not types even if they don't have any data
# that's not already in the type.
# This is also to be consistent with `StarAlgebras.MultiplicativeStructure`
(ordering::AbstractMonomialOrdering)(i, j) = cmp(ordering, i, j) < 0

"""
    struct LexOrder <: AbstractMonomialOrdering end

Lexicographic (Lex for short) Order often abbreviated as *lex* order as defined in [CLO13, Definition 2.2.3, p. 56]

The [`Graded`](@ref) version is often abbreviated as *grlex* order and is defined in [CLO13, Definition 2.2.5, p. 58]

[CLO13] Cox, D., Little, J., & OShea, D.
*Ideals, varieties, and algorithms: an introduction to computational algebraic geometry and commutative algebra*.
Springer Science & Business Media, **2013**.
"""
struct LexOrder <: AbstractMonomialOrdering end

const _TupleOrVector = Union{Tuple,AbstractVector}

function Base.cmp(::LexOrder, exp1::_TupleOrVector, exp2::_TupleOrVector)
    return cmp(exp1, exp2)
end

"""
    struct InverseLexOrder <: AbstractMonomialOrdering end

Inverse Lex Order defined in [CLO13, Exercise 2.2.6, p. 61] where it is abbreviated as *invlex*.
It corresponds to [`LexOrder`](@ref) but with the variables in reverse order.

The [`Graded`](@ref) version can be abbreviated as *grinvlex* order.
It is defined in [BDD13, Definition 2.1] where it is called *Graded xel order*.

[CLO13] Cox, D., Little, J., & OShea, D.
*Ideals, varieties, and algorithms: an introduction to computational algebraic geometry and commutative algebra*.
Springer Science & Business Media, **2013**.
[BDD13] Batselier, K., Dreesen, P., & De Moor, B.
*The geometry of multivariate polynomial division and elimination*.
SIAM Journal on Matrix Analysis and Applications, 34(1), 102-125, *2013*.
"""
struct InverseLexOrder <: AbstractMonomialOrdering end

# We can't use `Iterators.Reverse` because it's not an `AbstractVector`
# so not `cmp` methods is defined for it.
_rev(v::AbstractVector) = view(v, lastindex(v):-1:firstindex(v))
_rev(t::Tuple) = reverse(t)
function Base.cmp(::InverseLexOrder, exp1::_TupleOrVector, exp2::_TupleOrVector)
    return cmp(_rev(exp1), _rev(exp2))
end

"""
    struct Graded{O<:AbstractMonomialOrdering} <: AbstractMonomialOrdering
        same_degree_ordering::O
    end

Monomial ordering defined by:
* `degree(a) == degree(b)` then the ordering is determined by `same_degree_ordering`,
* otherwise, it is the ordering between the integers `degree(a)` and `degree(b)`.
"""
struct Graded{O<:AbstractMonomialOrdering} <: AbstractMonomialOrdering
    same_degree_ordering::O
end
Graded{O}() where {O<:AbstractMonomialOrdering} = Graded{O}(O())

function Base.cmp(ordering::Graded, a::_TupleOrVector, b::_TupleOrVector)
    deg_a = sum(a)
    deg_b = sum(b)
    if deg_a == deg_b
        return cmp(ordering.same_degree_ordering, a, b)
    else
        return deg_a - deg_b
    end
end

"""
    struct Reverse{O<:AbstractMonomialOrdering} <: AbstractMonomialOrdering
        reverse_order::O
    end

Monomial ordering defined by
`cmp(o::Reverse, a, b) where {O} = cmp(o.reverse_order, b, a)`.

Reverse Lex Order defined in [CLO13, Exercise 2.2.9, p. 61] where it is abbreviated as *rinvlex*.
can be obtained as `Reverse{InverseLexOrder}`.

The Graded Reverse Lex Order often abbreviated as *grevlex* order defined in [CLO13, Definition 2.2.6, p. 58]
can be obtained as `Graded{Reverse{InverseLexOrder}}`.

[CLO13] Cox, D., Little, J., & OShea, D.
*Ideals, varieties, and algorithms: an introduction to computational algebraic geometry and commutative algebra*.
Springer Science & Business Media, **2013**.
"""
struct Reverse{O<:AbstractMonomialOrdering} <: AbstractMonomialOrdering
    reverse_ordering::O
end
Reverse{O}() where {O<:AbstractMonomialOrdering} = Reverse{O}(O())

function Base.cmp(ordering::Reverse, a::_TupleOrVector, b::_TupleOrVector)
    return cmp(ordering.reverse_ordering, b, a)
end

#TODO(breaking) Return an instance, not a type
"""
    ordering(p::AbstractPolynomialLike)::Type{<:AbstractMonomialOrdering}

Returns the [`AbstractMonomialOrdering`](@ref) type to be used to compare
exponent vectors for the monomials of `p`.
"""
function ordering end

ordering(::Type{<:AbstractMonomial}) = Graded{LexOrder}
ordering(::Type{P}) where {P} = ordering(monomial_type(P))
ordering(p::AbstractPolynomialLike) = ordering(typeof(p))
# Useful for instance to ask ordering given the list
# of variables
ordering(::AbstractVector{T}) where {T} = ordering(T)
ordering(t::Tuple) = ordering(first(t))

# We reverse the order of comparisons here so that the result
# of x < y is equal to the result of Monomial(x) < Monomial(y)
# Without `Base.@pure`, TypedPolynomials allocates on Julia v1.6
# with `promote(x * y, x)`
Base.@pure function Base.cmp(v1::AbstractVariable, v2::AbstractVariable)
    return -cmp(name(v1), name(v2))
end

function Base.cmp(m1::AbstractMonomial, m2::AbstractMonomial)
    s1, s2 = promote_variables(m1, m2)
    return cmp(ordering(m1)(), exponents(s1), exponents(s2))
end

function compare(
    m1::AbstractMonomial,
    m2::AbstractMonomial,
    ::Type{O},
) where {O<:AbstractMonomialOrdering}
    s1, s2 = promote_variables(m1, m2)
    return cmp(O(), exponents(s1), exponents(s2))
end

# Implement this to make coefficients be compared with terms.
function _cmp_coefficient(a::Real, b::Real)
    return cmp(a, b)
end
function _cmp_coefficient(a::Number, b::Number)
    return cmp(abs(a), abs(b))
end
# By default, coefficients are not comparable so `a` is not strictly
# less than `b`, they are considered sort of equal.
_cmp_coefficient(a, b) = 0

function Base.cmp(t1::AbstractTermLike, t2::AbstractTermLike)
    Δ = cmp(monomial(t1), monomial(t2))
    if iszero(Δ)
        return _cmp_coefficient(coefficient(t1), coefficient(t2))
    end
    return Δ
end

Base.isless(t1::AbstractTermLike, t2::AbstractTermLike) = compare(t1, t2) < 0

"""
    struct ExponentsIterator{M}(
        object;
        mindegree::Int = 0,
        maxdegree::Union{Nothing,Int} = nothing,
        inline::Bool = false,
    )

An iterator for generating monomial exponents for monomial
ordering `M`. The type of the vector of exponents is the type of
`object` and is length (i.e., the number of variables) is `length(object)`.

Note that `object` does not have to be zero, it just needs to implement
`copy` and `setindex!` methods (except for `Tuple` which we handle with a
special case).

See also [`monomials`](@ref).

### Examples

The following example shows how to generate all exponents of
monomials of 2 variables up to degree 2.
```jldoctest
julia> collect(ExponentsIterator{Graded{LexOrder}}((0, 0), maxdegree = 2))
6-element Vector{Tuple{Int64, Int64}}:
 (0, 0)
 (0, 1)
 (1, 0)
 (0, 2)
 (1, 1)
 (2, 0)
```
Note that you can easily generate the tuple of exponents
of arbitrary length using `ntuple` as follows:
```jldoctest
julia> collect(ExponentsIterator{Graded{LexOrder}}(ntuple(zero, 3), mindegree = 2, maxdegree = 2))
6-element Vector{Tuple{Int64, Int64, Int64}}:
 (0, 0, 2)
 (0, 1, 1)
 (0, 2, 0)
 (1, 0, 1)
 (1, 1, 0)
 (2, 0, 0)
```
You can also change the monomial ordering and use `Vector` instead of `Tuple` as follows:
```jldoctest
julia> collect(ExponentsIterator{LexOrder}(zeros(Int, 2), mindegree = 2, maxdegree = 3))
7-element Vector{Vector{Int64}}:
 [0, 2]
 [0, 3]
 [1, 1]
 [1, 2]
 [2, 0]
 [2, 1]
 [3, 0]
```
"""
struct ExponentsIterator{M,D<:Union{Nothing,Int},O}
    object::O # Used to get number of variables and get new zero elements
    mindegree::Int
    maxdegree::D
    inline::Bool

    function ExponentsIterator{M}(
        object;
        mindegree::Int = 0,
        maxdegree::Union{Nothing,Int} = nothing,
        inline::Bool = false,
    ) where {M}
        if mindegree < 0
            throw(
                ArgumentError(
                    "The `mindegree` of `ExponentsIterator` cannot be negative.",
                ),
            )
        end
        if M <: Reverse
            throw(
                ArgumentError(
                    "Ordering `$M` is not a valid ordering, use `Graded{$M}` instead.",
                ),
            )
        end
        return new{M,typeof(maxdegree),typeof(object)}(
            object,
            mindegree,
            maxdegree,
            inline,
        )
    end
end

Base.eltype(::Type{ExponentsIterator{M,D,O}}) where {M,D,O} = O
# `IteratorSize` may actually be finite if the list of variables is empty.
# The same issue happens for `Iterators.Cycle`:
# https://github.com/JuliaLang/julia/pull/54187
# With `Base.IsInfinite`, some tests fail but with `Base.SizeUnknown`, all
# tests pass.
function Base.IteratorSize(::Type{<:ExponentsIterator{M,Nothing}}) where {M}
    return Base.SizeUnknown()
end
function Base.IteratorSize(::Type{<:ExponentsIterator{M,Int}}) where {M}
    return Base.HasLength()
end

function _length(it::ExponentsIterator, maxdegree)
    if maxdegree < it.mindegree
        return 0
    end
    len = binomial(nvariables(it) + maxdegree, nvariables(it))
    if it.mindegree > 0
        len -= binomial(nvariables(it) + it.mindegree - 1, nvariables(it))
    end
    return len
end

function Base.length(it::ExponentsIterator{M,Int}) where {M}
    return _length(it, it.maxdegree)
end

function Base.length(it::ExponentsIterator{M,Nothing}) where {M}
    if isempty(it.object)
        return _length(it, it.mindegree)
    else
        error("The iterator is infinity because `maxdegree` is `nothing`.")
    end
end

nvariables(it::ExponentsIterator) = length(it.object)

_last_lex_index(n, ::Type{LexOrder}) = n
_prev_lex_index(i, ::Type{LexOrder}) = i - 1
_not_first_indices(n, ::Type{LexOrder}) = n:-1:2
_last_lex_index(_, ::Type{InverseLexOrder}) = 1
_prev_lex_index(i, ::Type{InverseLexOrder}) = i + 1
_not_first_indices(n, ::Type{InverseLexOrder}) = 1:(n-1)
_last_lex_index(n, ::Type{Reverse{LexOrder}}) = 1
_prev_lex_index(i, ::Type{Reverse{LexOrder}}) = i + 1
_not_first_indices(n, ::Type{Reverse{LexOrder}}) = (n-1):-1:1
_last_lex_index(n, ::Type{Reverse{InverseLexOrder}}) = n
_prev_lex_index(i, ::Type{Reverse{InverseLexOrder}}) = i - 1
_not_first_indices(n, ::Type{Reverse{InverseLexOrder}}) = 2:n
_last_lex_index(n, ::Type{Graded{M}}) where {M} = _last_lex_index(n, M)
_prev_lex_index(i, ::Type{Graded{M}}) where {M} = _prev_lex_index(i, M)
_not_first_indices(n, ::Type{Graded{M}}) where {M} = _not_first_indices(n, M)

_increase_degree(it::ExponentsIterator{<:Graded,Nothing}, _) = false
_increase_degree(it::ExponentsIterator{<:Graded,Int}, _) = false
_increase_degree(it::ExponentsIterator{M,Nothing}, _) where {M} = true
function _increase_degree(it::ExponentsIterator{M,Int}, deg) where {M}
    return deg < it.maxdegree
end

# We just changed the degree by removing `Δ`,
# In graded ordering, we just add `Δ` to maintain the same degree
_adjust_degree(::ExponentsIterator{<:Graded}, _, Δ) = Δ
# Otherwise, we just need the degree to stay above `it.mindegree`,
# so we need to add `it.mindegree - deg`
_adjust_degree(it::ExponentsIterator, deg, _) = max(0, it.mindegree - deg)

# Same as `BangBang.setindex!!`
_setindex!(x, v, i) = Base.setindex!(x, v, i)
_setindex!(x::Tuple, v, i) = Base.setindex(x, v, i)
_increment!(x, Δ, i) = _setindex!(x, x[i] + Δ, i)

_zero(x) = zero(x)
_zero(x::Tuple) = zero.(x)

_zero!(x) = fill!(x, 0)
_zero!(x::Tuple) = _zero(x)

_copy(x) = copy(x)
_copy(x::Tuple) = x

function __iterate!(
    it::ExponentsIterator{Graded{Reverse{M}}},
    z,
    i,
    deg,
) where {M}
    z = _increment!(z, -1, i)
    k = _last_lex_index(nvariables(it), M)
    Δ = z[k] + 1
    z = _setindex!(z, 0, k)
    z = _setindex!(z, Δ, _prev_lex_index(i, Reverse{M}))
    return z, deg
end

function __iterate!(it::ExponentsIterator{M}, z, i, deg) where {M}
    Δ = z[i] - 1
    z = _setindex!(z, 0, i)
    deg -= Δ
    Δ = _adjust_degree(it, deg, Δ)
    deg += Δ
    z = _setindex!(z, Δ, _last_lex_index(nvariables(it), M))
    z = _increment!(z, 1, _prev_lex_index(i, M))
    return z, deg
end

function _iterate!(it::ExponentsIterator{M}, z, deg) where {M}
    if _increase_degree(it, deg)
        z = _increment!(z, 1, _last_lex_index(nvariables(it), M))
        return z, deg + 1
    end
    I = _not_first_indices(nvariables(it), M)
    i = findfirst(i -> !iszero(z[i]), I)
    if isnothing(i)
        if !isnothing(it.maxdegree) && deg == it.maxdegree
            return
        end
        z = _zero!(z)
        z = _setindex!(z, deg + 1, _last_lex_index(nvariables(it), M))
        return z, deg + 1
    end
    return __iterate!(it, z, I[i], deg)
end

function Base.iterate(it::ExponentsIterator{M}) where {M}
    z = _zero(it.object)
    if !isnothing(it.maxdegree) && it.maxdegree < it.mindegree
        return
    end
    if it.mindegree > 0
        if nvariables(it) == 0
            return
        end
        z = _setindex!(z, it.mindegree, _last_lex_index(nvariables(it), M))
    end
    return z, (z, it.mindegree)
end

function Base.iterate(it::ExponentsIterator, state)
    if nvariables(it) == 0
        return # There cannot be more than 1 element
    end
    z, deg = state
    if !it.inline
        z = _copy(z)
    end
    state = _iterate!(it, z, deg)
    if isnothing(state)
        return
    end
    return state[1], state
end

# TODO Backward compat, remove the following in next breaking release
"""
    compare(a, b, order::Type{<:AbstractMonomialOrdering})

Returns a negative number if `a < b`, a positive number if `a > b` and zero if `a == b`.
The comparison is done according to `order`.

**Warning** This is deprecated, use `cmp(order(), a, b)` instead.
"""
function compare end

function compare(t1::AbstractTermLike, t2::AbstractTermLike)
    return compare(t1, t2, ordering(t1))
end

function compare(
    e1::_TupleOrVector,
    e2::_TupleOrVector,
    ::Type{O},
) where {O<:AbstractMonomialOrdering}
    return cmp(O(), e1, e2)
end

function compare(
    t1::AbstractTermLike,
    t2::AbstractTermLike,
    ::Type{O},
) where {O<:AbstractMonomialOrdering}
    Δ = compare(monomial(t1), monomial(t2), O)
    if iszero(Δ)
        return _cmp_coefficient(coefficient(t1), coefficient(t2))
    end
    return Δ
end

function compare(
    a::AbstractMonomial,
    b::AbstractMonomial,
    ::Type{Graded{O}},
) where {O}
    deg_a = degree(a)
    deg_b = degree(b)
    if deg_a == deg_b
        return compare(a, b, O)
    else
        return deg_a - deg_b
    end
end

function compare(
    a::AbstractMonomial,
    b::AbstractMonomial,
    ::Type{Reverse{O}},
) where {O}
    return compare(b, a, O)
end

Base.isless(v1::AbstractVariable, v2::AbstractVariable) = cmp(v1, v2) < 0

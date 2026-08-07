"""
    promote_variables(p::AbstractPolynomialLike, q::AbstractPolynomialLike)

Return two polynomials over the same variables.
"""
function promote_variables end

# Simplified promote rules: Julia's `promote` tries both orderings automatically,
# so we only define one direction per pair.

# MonomialLike: promote_rule delegates to monomial_type promotion
Base.promote_rule(::Type{M}, ::Type{M}) where {M<:AbstractMonomialLike} = M
function Base.promote_rule(
    M1::Type{<:AbstractMonomialLike},
    M2::Type{<:AbstractMonomialLike},
)
    return promote_type(monomial_type(M1), monomial_type(M2))
end

# SA.Term: promote coefficient and monomial types
function Base.promote_rule(
    TS::Type{<:SA.Term{S}},
    TT::Type{<:SA.Term{T}},
) where {S,T}
    U = promote_type(S, T)
    M = promote_type(monomial_type(TS), monomial_type(TT))
    return term_type(M, U)
end

# MonomialLike + Term
function Base.promote_rule(
    TS::Type{<:AbstractMonomialLike},
    TT::Type{<:SA.Term{T}},
) where {T}
    U = promote_type(Int, T)
    M = promote_type(monomial_type(TS), monomial_type(TT))
    return term_type(M, U)
end

# PolynomialLike (_APL): promote via term_type/polynomial_type
Base.promote_rule(::Type{PT}, ::Type{PT}) where {PT<:_APL} = PT
function Base.promote_rule(PS::Type{<:_APL}, PT::Type{<:_APL})
    TS = try
        term_type(PS)
    catch
        return Union{}
    end
    TT = try
        term_type(PT)
    catch
        return Union{}
    end
    return polynomial_type(promote_type(TS, TT))
end

# Constant (Number) × Polynomial
function promote_rule_constant(::Type{S}, PT::Type{<:_APL{T}}) where {S<:Number,T}
    try
        return polynomial_type(PT, promote_type(S, T))
    catch
        return Any
    end
end
function Base.promote_rule(::Type{PT}, ::Type{T}) where {T<:Number,PT<:_APL}
    return promote_rule_constant(T, PT)
end

# Rational
function promote_rule_constant(
    ::Type{T},
    ::Type{RationalPoly{NT,DT}},
) where {T,NT,DT}
    return RationalPoly{promote_type(T, NT),promote_type(DT, term_type(DT))}
end

function Base.promote_rule(::Type{RT}, ::Type{T}) where {T,RT<:RationalPoly}
    return promote_rule_constant(T, RT)
end

function promote_rule_rational(
    ::Type{PT},
    ::Type{RationalPoly{S,T}},
) where {PT<:_APL,S,T}
    return RationalPoly{promote_type(PT, S),promote_type(T, term_type(T))}
end
function promote_rule_rational(
    ::Type{RationalPoly{S,T}},
    ::Type{RationalPoly{U,V}},
) where {S,T,U,V}
    return RationalPoly{promote_type(S, U),promote_type(T, V)}
end

function Base.promote_rule(
    ::Type{RS},
    ::Type{RT},
) where {RS<:RationalPoly,RT<:RationalPoly}
    return promote_rule_rational(RS, RT)
end
function Base.promote_rule(
    ::Type{PT},
    ::Type{RT},
) where {PT<:_APL,RT<:RationalPoly}
    return promote_rule_rational(PT, RT)
end

# MutableArithmetics
function MA.promote_operation(
    op::Union{typeof(+),typeof(-)},
    PT::Type{<:_APL{S}},
    QT::Type{<:_APL{T}},
) where {S,T}
    U = MA.promote_operation(op, S, T)
    return polynomial_type(
        promote_type(monomial_type(PT), monomial_type(QT)),
        U,
    )
end
function MA.promote_operation(
    ::typeof(*),
    MT1::Type{<:AbstractMonomialLike},
    MT2::Type{<:AbstractMonomialLike},
)
    return promote_type(monomial_type(MT1), monomial_type(MT2))
end
function MA.promote_operation(
    ::typeof(*),
    TT::Type{<:SA.Term{S}},
    ST::Type{<:SA.Term{T}},
) where {S,T}
    UT = MA.promote_operation(*, monomial_type(TT), monomial_type(ST))
    U = MA.promote_operation(*, S, T)
    return promote_operation_left_constant(*, U, UT)
end
function MA.promote_operation(
    ::typeof(*),
    TT::Type{<:AbstractMonomialLike},
    ST::Type{<:SA.Term{T}},
) where {T}
    UT = MA.promote_operation(*, monomial_type(TT), monomial_type(ST))
    U = MA.promote_operation(*, Int, T)
    return promote_operation_left_constant(*, U, UT)
end
function MA.promote_operation(
    ::typeof(*),
    PT::Type{<:_APL{S}},
    QT::Type{<:_APL{T}},
) where {S,T}
    UP = MA.promote_operation(*, monomial_type(PT), monomial_type(QT))
    U = MA.promote_sum_mul(S, T)
    return polynomial_type(promote_operation_left_constant(*, U, UP))
end

function promote_operation_left_constant(
    ::typeof(*),
    ::Type{T},
    ::Type{M},
) where {T,M<:AbstractMonomialLike}
    return term_type(M, T)
end

function promote_operation_right_constant(
    ::typeof(*),
    ::Type{M},
    ::Type{T},
) where {T,M<:AbstractMonomialLike}
    return term_type(M, T)
end

function promote_operation_left_constant(
    ::typeof(*),
    ::Type{T},
    ::Type{P},
) where {T,U,P<:_APL{U}}
    return similar_type(P, MA.promote_operation(*, T, U))
end

function promote_operation_right_constant(
    ::typeof(*),
    ::Type{P},
    ::Type{T},
) where {T,U,P<:_APL{U}}
    return similar_type(P, MA.promote_operation(*, U, T))
end

function MA.promote_operation(
    ::typeof(*),
    ::Type{T},
    ::Type{P},
) where {T<:Number,P<:_APL}
    return promote_operation_left_constant(*, T, P)
end

function MA.promote_operation(
    ::typeof(*),
    ::Type{P},
    ::Type{T},
) where {T<:Number,P<:_APL}
    return promote_operation_right_constant(*, P, T)
end

function MA.promote_operation(
    ::typeof(*),
    ::Type{P},
    ::Type{RationalPoly{NT,DT}},
) where {P<:_APL,NT,DT}
    return RationalPoly{MA.promote_operation(*, P, NT),DT}
end

function MA.promote_operation(
    ::typeof(*),
    ::Type{RationalPoly{NT,DT}},
    ::Type{P},
) where {P<:_APL,NT,DT}
    return RationalPoly{MA.promote_operation(*, NT, P),DT}
end

function MA.promote_operation(
    ::typeof(*),
    ::Type{RationalPoly{NS,DS}},
    ::Type{RationalPoly{NT,DT}},
) where {NS,DS,NT,DT}
    return RationalPoly{
        MA.promote_operation(*, NS, NT),
        MA.promote_operation(*, DS, DT),
    }
end

function _search_sorted_first(haystack::AbstractVector, needle; kws...)
    return searchsortedfirst(haystack, needle; kws...)
end
function _search_sorted_first(haystack::Tuple, needle; kws...)
    return findfirst(isequal(needle), haystack)
end

_idx(needle, haystack) = _search_sorted_first(haystack, needle, rev = true)

struct ExponentMap{I,L} <: Function
    indices::I
    length::L
end

function (map::ExponentMap{Vector{Int}})(exp::Vector{Int})
    new_exp = zeros(Int, map.length)
    for (i, e) in zip(map.indices, exp)
        new_exp[i] = e
    end
    return new_exp
end

function (map::ExponentMap{NTuple{N,Int}})(exp::NTuple{N,Int}) where {N}
    return ntuple(map.length::Val) do i
        j = findfirst(isequal(i), map.indices)
        if isnothing(j)
            return 0
        else
            return exp[j]
        end
    end
end

_length(x::AbstractVector) = length(x)
_length(::NTuple{N,Any}) where {N} = Val(N)

function _map(needles, haystack)
    if length(needles) == length(haystack)
        return nothing
    end
    return ExponentMap(
        map(Base.Fix2(_idx, haystack), needles),
        _length(haystack),
    )
end

"""
    promote_variables_with_maps(a, b)

Given two sorted variable collections `a` and `b`, return
`((all_vars, map_a), (all_vars, map_b))` where `all_vars` is the merged
sorted variable set and `map_a` (resp. `map_b`) is an `ExponentMap`
that maps exponents from `a` (resp. `b`) to exponents in `all_vars`,
or `nothing` if no mapping is needed.
"""
function promote_variables_with_maps(a, b)
    if a == b
        return (a, nothing), (b, nothing)
    end
    all_vars = SA.merge_sorted(
        a,
        b;
        lt = isless,
        combine = SA.first_of,
        filter = _ -> true,
        rev = true,
    )
    return (all_vars, _map(a, all_vars)), (all_vars, _map(b, all_vars))
end

function SA.promote_with_map(t::AbstractTerm, all_vars, map::ExponentMap)
    mono, _ = SA.promote_with_map(monomial(t), all_vars, map)
    return term(coefficient(t), mono), map
end

function SA.promote_with_map(p::AbstractPolynomial, all_vars, map::ExponentMap)
    new_terms = [
        term(
            coefficient(t),
            first(SA.promote_with_map(monomial(t), all_vars, map)),
        ) for t in terms(p)
    ]
    return polynomial(new_terms, SortedUniqState()), map
end

function SA.promote_bases_with_maps(p::_APL, q::_APL)
    _p, _q = promote_variables_with_maps(variables(p), variables(q))
    return SA.maybe_promote(p, _p...), SA.maybe_promote(q, _q...)
end

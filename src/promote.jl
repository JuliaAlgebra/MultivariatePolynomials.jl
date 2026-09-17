function promote_variables end

Base.promote_rule(::Type{M}, ::Type{M}) where {M<:AbstractMonomialLike} = M
function Base.promote_rule(
    ::Type{M},
    ::Type{N},
) where {M<:AbstractMonomialLike,N<:AbstractMonomialLike}
    return promote_type(monomial_type(M), monomial_type(N))
end
function MA.promote_operation(
    ::typeof(*),
    ::Type{M},
    ::Type{N},
) where {M<:AbstractMonomialLike,N<:AbstractMonomialLike}
    return promote_type(monomial_type(M), monomial_type(N))
end
function MA.promote_operation(
    op::Union{typeof(+),typeof(-)},
    ::Type{M},
    ::Type{N},
) where {M<:AbstractMonomialLike,N<:AbstractMonomialLike}
    return polynomial_type(MA.promote_operation(*, M, N))
end
function MA.promote_operation(
    ::typeof(*),
    ::Type{T},
    ::Type{M},
) where {T<:Number,M<:AbstractMonomialLike}
    return term_type(M, T)
end
function MA.promote_operation(
    ::typeof(*),
    ::Type{M},
    ::Type{T},
) where {M<:AbstractMonomialLike,T<:Number}
    return term_type(M, T)
end

function MA.promote_operation(
    op::Union{typeof(+),typeof(-),typeof(*)},
    ::Type{T},
    ::Type{S},
) where {T<:AbstractTerm,S<:AbstractTerm}
    return MA.promote_operation(op, polynomial_type(T), polynomial_type(S))
end

# Rational
function promote_rule_constant(
    ::Type{T},
    ::Type{RationalPoly{NT,DT}},
) where {T<:Number,NT,DT}
    return RationalPoly{promote_type(T, NT),promote_type(DT, term_type(DT))}
end

function Base.promote_rule(
    ::Type{RT},
    ::Type{T},
) where {T<:Number,RT<:RationalPoly}
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

function SA.promote_with_map(p::Polynomial{B}, vars, map::ExponentMap) where {B}
    return Polynomial(Variables{B}(vars), map(exponents(p))), map
end
function SA.promote_bases_with_maps(
    p::Polynomial{B},
    q::Polynomial{B},
) where {B}
    a, b = promote_variables_with_maps(variables(p), variables(q))
    return SA.maybe_promote(p, a...), SA.maybe_promote(q, b...)
end
function promote_variables(p::AbstractMonomial, q::AbstractMonomial)
    return SA.promote_bases(p, q)
end

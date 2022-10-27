# MonomialLike
Base.promote_rule(::Type{M}, ::Type{M}) where {M<:AbstractMonomialLike} = M
function Base.promote_rule(M1::Type{<:AbstractMonomialLike},
                           M2::Type{<:AbstractMonomialLike})
    return promote_type(monomialtype(M1), monomialtype(M2))
end

# TermLike
Base.promote_rule(::Type{T}, ::Type{T}) where {T<:AbstractTermLike} = T
function Base.promote_rule(TS::Type{<:AbstractTermLike{S}}, TT::Type{<:AbstractTermLike{T}}) where {S, T}
    U = promote_type(S, T)
    M = promote_type(monomialtype(TS), monomialtype(TT))
    return termtype(M, U)
end
function promote_rule_constant(::Type{S}, TT::Type{<:AbstractTermLike{T}}) where {S, T}
    return termtype(TT, promote_type(S, T))
end


# PolynomialLike
Base.promote_rule(::Type{PT}, ::Type{PT}) where {PT<:APL} = PT
function Base.promote_rule(PS::Type{<:APL}, PT::Type{<:APL})
    return polynomialtype(promote_type(termtype(PS), termtype(PT)))
end

function promote_rule_constant(::Type{S}, PT::Type{<:APL{T}}) where {S, T}
    return polynomialtype(PT, promote_type(S, T))
end
Base.promote_rule(::Type{T}, ::Type{PT}) where {T, PT<:APL} = promote_rule_constant(T, PT)
# Resolve method ambiguity with Base:
Base.promote_rule(::Type{Any}, ::Type{<:APL}) = Any

# We don't have any information on the MultivariatePolynomials implementation,
# so we won't be able to convert the constant to `APL`.
promote_rule_constant(::Type, PT::Type{AbstractMonomialLike}) = Any
promote_rule_constant(::Type, PT::Type{AbstractTermLike{T}}) where {T} = Any
promote_rule_constant(::Type, PT::Type{AbstractTermLike}) = Any
promote_rule_constant(::Type, PT::Type{APL{T}}) where {T} = Any
promote_rule_constant(::Type, PT::Type{APL}) = Any

# AbstractMonomialLike{T}
Base.promote_rule(::Type{AbstractMonomialLike}, ::Type{<:AbstractMonomialLike}) = AbstractMonomialLike
Base.promote_rule(::Type{<:AbstractMonomialLike}, ::Type{AbstractMonomialLike}) = AbstractMonomialLike
Base.promote_rule(::Type{AbstractMonomialLike}, ::Type{<:AbstractTermLike{T}}) where {T} = _atl(Int, T)
Base.promote_rule(::Type{<:AbstractTermLike{T}}, ::Type{AbstractMonomialLike}) where {T} = _atl(Int, T)
Base.promote_rule(::Type{AbstractMonomialLike}, ::Type{AbstractTermLike{T}}) where {T} = _atl(Int, T)
Base.promote_rule(::Type{AbstractTermLike{T}}, ::Type{AbstractMonomialLike}) where {T} = _atl(Int, T)
Base.promote_rule(::Type{AbstractMonomialLike}, ::Type{<:APL{T}}) where {T} = _apl(Int, T)
Base.promote_rule(::Type{<:APL{T}}, ::Type{AbstractMonomialLike}) where {T} = _apl(Int, T)
Base.promote_rule(::Type{AbstractMonomialLike}, ::Type{APL{T}}) where {T} = _apl(Int, T)
Base.promote_rule(::Type{APL{T}}, ::Type{AbstractMonomialLike}) where {T} = _apl(Int, T)

# AbstractTermLike{T}
_atl(::Type{T}, ::Type{T}) where {T} = AbstractTermLike{T}
_atl(::Type, ::Type) = AbstractTermLike
__atl(::Type{T}, ::Type{<:AbstractTermLike{S}}) where {S,T} = _atl(T, S)
__atl(::Type{T}, ::Type{<:APL{S}}) where {S,T} = _apl(T, S)
Base.promote_rule(::Type{AbstractTermLike{T}}, P::Type{<:AbstractTermLike{S}}) where {S,T} = _atl(T, S)
Base.promote_rule(P::Type{<:AbstractTermLike{S}}, ::Type{AbstractTermLike{T}}) where {S,T} = _atl(T, S)
Base.promote_rule(::Type{AbstractTermLike{T}}, P::Type{<:APL{S}}) where {S,T} = _apl(T, S)
Base.promote_rule(P::Type{<:APL{S}}, ::Type{AbstractTermLike{T}}) where {S,T} = _apl(T, S)
Base.promote_rule(::Type{AbstractTermLike{T}}, P::Type{APL{S}}) where {S,T} = _apl(T, S)
Base.promote_rule(P::Type{APL{S}}, ::Type{AbstractTermLike{T}}) where {S,T} = _apl(T, S)

# AbstractTermLike
Base.promote_rule(::Type{AbstractTermLike}, ::Type{<:AbstractTermLike}) = AbstractTermLike
Base.promote_rule(::Type{<:AbstractTermLike}, ::Type{AbstractTermLike}) = AbstractTermLike
Base.promote_rule(::Type{AbstractTermLike}, ::Type{<:APL}) = APL
Base.promote_rule(::Type{<:APL}, ::Type{AbstractTermLike}) = APL
Base.promote_rule(::Type{AbstractTermLike}, ::Type{APL}) = APL
Base.promote_rule(::Type{APL}, ::Type{AbstractTermLike}) = APL

# APL{T}
_apl(::Type{T}, ::Type{T}) where {T} = APL{T}
_apl(::Type, ::Type) = APL
Base.promote_rule(::Type{APL{T}}, ::Type{<:APL{S}}) where {S,T} = _apl(S, T)
Base.promote_rule(::Type{<:APL{S}}, ::Type{APL{T}}) where {S,T} = _apl(S, T)

# APL
Base.promote_rule(::Type{APL}, ::Type{<:APL}) = APL
Base.promote_rule(::Type{<:APL}, ::Type{APL}) = APL

# Rational
promote_rule_constant(::Type{T}, ::Type{RationalPoly{NT, DT}}) where {T, NT, DT} = RationalPoly{promote_type(T, NT), promote_type(DT, termtype(DT))}

Base.promote_rule(::Type{T}, ::Type{RT}) where {T, RT<:RationalPoly} = promote_rule_constant(T, RT)
Base.promote_rule(::Type{RT}, ::Type{T}) where {T, RT<:RationalPoly} = promote_rule_constant(T, RT)

promote_rule_rational(::Type{PT}, ::Type{RationalPoly{S, T}}) where {PT<:APL, S, T} = RationalPoly{promote_type(PT, S), promote_type(T, termtype(T))}
promote_rule_rational(::Type{RationalPoly{S, T}}, ::Type{RationalPoly{U, V}}) where {S, T, U, V} = RationalPoly{promote_type(S, U), promote_type(T, V)}

Base.promote_rule(::Type{RS}, ::Type{RT}) where {RS<:RationalPoly, RT<:RationalPoly} = promote_rule_rational(RS, RT)
Base.promote_rule(::Type{PT}, ::Type{RT}) where {PT<:APL, RT<:RationalPoly} = promote_rule_rational(PT, RT)
Base.promote_rule(::Type{RT}, ::Type{PT}) where {PT<:APL, RT<:RationalPoly} = promote_rule_rational(PT, RT)

# MutableArithmetics
function MA.promote_operation(
    op::Union{typeof(+), typeof(-)}, PT::Type{<:APL{S}},
    QT::Type{<:APL{T}}) where {S, T}

    U = MA.promote_operation(op, S, T)
    return polynomialtype(promote_type(monomialtype(PT), monomialtype(QT)), U)
end
function MA.promote_operation(::typeof(*), MT1::Type{<:AbstractMonomialLike},
    MT2::Type{<:AbstractMonomialLike})
    return promote_type(monomialtype(MT1), monomialtype(MT2))
end
function MA.promote_operation(::typeof(*), TT::Type{<:AbstractTermLike{S}},
                              ST::Type{<:AbstractTermLike{T}}) where {S, T}
    U = MA.promote_operation(*, S, T)
    return termtype(promote_type(monomialtype(TT), monomialtype(ST)), U)
end
function MA.promote_operation(::typeof(*), PT::Type{<:APL{S}}, QT::Type{<:APL{T}}) where {S, T}
    ST = MA.promote_operation(*, S, T)
    U = MA.promote_operation(+, ST, ST)
    return polynomialtype(promote_type(monomialtype(PT), monomialtype(QT)), U)
end
function MA.promote_operation(::typeof(*), ::Type{T}, ::Type{P}) where {T, U, P<:APL{U}}
    return changecoefficienttype(P, MA.promote_operation(*, T, U))
end
function MA.promote_operation(::typeof(*), ::Type{P}, ::Type{T}) where {T, U, P<:APL{U}}
    return changecoefficienttype(P, MA.promote_operation(*, U, T))
end

# MonomialLike
Base.promote_rule(::Type{M}, ::Type{M}) where {M<:AbstractMonomialLike} = M
function Base.promote_rule(M1::Type{<:AbstractMonomialLike},
                           M2::Type{<:AbstractMonomialLike})
    return promote_type(monomial_type(M1), monomial_type(M2))
end

# TermLike
Base.promote_rule(::Type{T}, ::Type{T}) where {T<:AbstractTermLike} = T
function Base.promote_rule(TS::Type{<:AbstractTermLike{S}}, TT::Type{<:AbstractTermLike{T}}) where {S, T}
    U = promote_type(S, T)
    M = promote_type(monomial_type(TS), monomial_type(TT))
    return term_type(M, U)
end
function promote_rule_constant(::Type{S}, TT::Type{<:AbstractTermLike{T}}) where {S, T}
    return term_type(TT, promote_type(S, T))
end


# PolynomialLike
Base.promote_rule(::Type{PT}, ::Type{PT}) where {PT<:APL} = PT
function Base.promote_rule(PS::Type{<:APL}, PT::Type{<:APL})
    return polynomial_type(promote_type(term_type(PS), term_type(PT)))
end

function promote_rule_constant(::Type{S}, PT::Type{<:APL{T}}) where {S, T}
    return polynomial_type(PT, promote_type(S, T))
end
Base.promote_rule(::Type{PT}, ::Type{T}) where {T, PT<:APL} = promote_rule_constant(T, PT)

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
promote_rule_constant(::Type{T}, ::Type{RationalPoly{NT, DT}}) where {T, NT, DT} = RationalPoly{promote_type(T, NT), promote_type(DT, term_type(DT))}

Base.promote_rule(::Type{RT}, ::Type{T}) where {T, RT<:RationalPoly} = promote_rule_constant(T, RT)

promote_rule_rational(::Type{PT}, ::Type{RationalPoly{S, T}}) where {PT<:APL, S, T} = RationalPoly{promote_type(PT, S), promote_type(T, term_type(T))}
promote_rule_rational(::Type{RationalPoly{S, T}}, ::Type{RationalPoly{U, V}}) where {S, T, U, V} = RationalPoly{promote_type(S, U), promote_type(T, V)}

Base.promote_rule(::Type{RS}, ::Type{RT}) where {RS<:RationalPoly, RT<:RationalPoly} = promote_rule_rational(RS, RT)
Base.promote_rule(::Type{PT}, ::Type{RT}) where {PT<:APL, RT<:RationalPoly} = promote_rule_rational(PT, RT)
Base.promote_rule(::Type{RT}, ::Type{PT}) where {PT<:APL, RT<:RationalPoly} = promote_rule_rational(PT, RT)

# MutableArithmetics
function MA.promote_operation(
    op::Union{typeof(+), typeof(-)}, PT::Type{<:APL{S}},
    QT::Type{<:APL{T}}) where {S, T}

    U = MA.promote_operation(op, S, T)
    return polynomial_type(promote_type(monomial_type(PT), monomial_type(QT)), U)
end
function MA.promote_operation(::typeof(*), MT1::Type{<:AbstractMonomialLike},
    MT2::Type{<:AbstractMonomialLike})
    return promote_type(monomial_type(MT1), monomial_type(MT2))
end
function MA.promote_operation(::typeof(*), TT::Type{<:AbstractTermLike{S}},
                              ST::Type{<:AbstractTermLike{T}}) where {S, T}
    U = MA.promote_operation(*, S, T)
    return term_type(promote_type(monomial_type(TT), monomial_type(ST)), U)
end
function MA.promote_operation(::typeof(*), PT::Type{<:APL{S}}, QT::Type{<:APL{T}}) where {S, T}
    ST = MA.promote_operation(*, S, T)
    U = MA.promote_operation(+, ST, ST)
    return polynomial_type(promote_type(monomial_type(PT), monomial_type(QT)), U)
end
function MA.promote_operation(::typeof(*), ::Type{T}, ::Type{P}) where {T, U, P<:APL{U}}
    return similar_type(P, MA.promote_operation(*, T, U))
end
function MA.promote_operation(::typeof(*), ::Type{P}, ::Type{T}) where {T, U, P<:APL{U}}
    return similar_type(P, MA.promote_operation(*, U, T))
end

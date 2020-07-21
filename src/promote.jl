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
Base.promote_rule(::Type{PT}, ::Type{T}) where {T, PT<:APL} = promote_rule_constant(T, PT)
Base.promote_rule(::Type{T}, ::Type{PT}) where {T, PT<:APL} = promote_rule_constant(T, PT)

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
function MA.promote_operation(::typeof(*), PT::Type{<:APL{S}},
                              QT::Type{<:APL{T}}) where {S, T}
    U = MA.promote_operation(MA.add_mul, S, T)
    return polynomialtype(promote_type(monomialtype(PT), monomialtype(QT)), U)
end

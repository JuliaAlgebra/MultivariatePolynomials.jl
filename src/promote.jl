Base.promote_rule(::Type{PT}, ::Type{PS}) where {PT<:APL, PS<:APL} = promote_type(polynomialtype(PT), polynomialtype(PS))
Base.promote_rule(::Type{PT}, ::Type{PT}) where {PT<:APL} = PT

promote_rule_constant(::Type{T}, ::Type{RationalPoly{NT, DT}}) where {T, NT, DT} = RationalPoly{promote_type(T, NT), promote_type(DT, termtype(DT))}

Base.promote_rule(::Type{PT}, ::Type{T}) where {T, PT<:APL} = promote_rule_constant(T, PT)
Base.promote_rule(::Type{T}, ::Type{PT}) where {T, PT<:APL} = promote_rule_constant(T, PT)
Base.promote_rule(::Type{T}, ::Type{RT}) where {T, RT<:RationalPoly} = promote_rule_constant(T, RT)
Base.promote_rule(::Type{RT}, ::Type{T}) where {T, RT<:RationalPoly} = promote_rule_constant(T, RT)

promote_rule_rational(::Type{PT}, ::Type{RationalPoly{S, T}}) where {PT<:APL, S, T} = RationalPoly{promote_type(PT, S), promote_type(T, termtype(T))}
promote_rule_rational(::Type{RationalPoly{S, T}}, ::Type{RationalPoly{U, V}}) where {S, T, U, V} = RationalPoly{promote_type(S, U), promote_type(T, V)}

Base.promote_rule(::Type{RS}, ::Type{RT}) where {RS<:RationalPoly, RT<:RationalPoly} = promote_rule_rational(RS, RT)
Base.promote_rule(::Type{PT}, ::Type{RT}) where {PT<:APL, RT<:RationalPoly} = promote_rule_rational(PT, RT)
Base.promote_rule(::Type{RT}, ::Type{PT}) where {PT<:APL, RT<:RationalPoly} = promote_rule_rational(PT, RT)

# Promotion with Term
function Base.promote_rule(ST::Type{<:AbstractTermLike{S}}, TT::Type{<:AbstractTerm{T}}) where {S, T}
    U = promote_type(S, T)
    UT = termtype(ST, U)
    if UT != termtype(TT, U)
        error("Cannot promote `$ST` and `$TT` to the same type.")
    end
    return UT
end

#promote_rule(::Type{Term{C, U}}, ::Type{RationalPoly{C, S, T}}) where {C, S, T, U} = RationalPoly{C, promote_type(U, S), T}
#promote_rule(::Type{RationalPoly{C, S, T}}, ::Type{Term{C, U}}) where {C, S, T, U} = RationalPoly{C, promote_type(U, S), T}
#promote_rule(::Type{Polynomial{C, U}}, ::Type{RationalPoly{C, S, T}}) where {C, S, T, U} = RationalPoly{C, promote_type(U, S), T}
#promote_rule(::Type{RationalPoly{C, S, T}}, ::Type{Polynomial{C, U}}) where {C, S, T, U} = RationalPoly{C, promote_type(U, S), T}

function MA.promote_operation(
    op::Union{typeof(+), typeof(-)}, PT::Type{<:APL{S}},
    QT::Type{<:APL{T}}) where {S, T}

    U = MA.promote_operation(op, S, T)
    return polynomialtype(promote_type(monomialtype(PT), monomialtype(QT)), U)
end
function MA.promote_operation(::typeof(*), MT1::Type{<:AbstractMonomialLike},
    MT2::Type{<:AbstractMonomialLike})
    return typeof(constantmonomial(MT1) * constantmonomial(MT2))
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

Base.promote_rule(::Type{PT}, ::Type{PS}) where {PT<:APL, PS<:APL} = promote_rule(polynomialtype(PT), polynomialtype(PS))
Base.promote_rule(::Type{PT}, ::Type{PT}) where {PT<:APL} = PT

promote_rule_constant(::Type{T}, ::Type{RationalPoly{NT, DT}}) where {T, NT, DT} = RationalPoly{promote_rule(T, NT), promote_rule(DT, termtype(DT))}

Base.promote_rule(::Type{PT}, ::Type{T}) where {T, PT<:APL} = promote_rule_constant(T, PT)
Base.promote_rule(::Type{T}, ::Type{PT}) where {T, PT<:APL} = promote_rule_constant(T, PT)
Base.promote_rule(::Type{T}, ::Type{RT}) where {T, RT<:RationalPoly} = promote_rule_constant(T, RT)
Base.promote_rule(::Type{RT}, ::Type{T}) where {T, RT<:RationalPoly} = promote_rule_constant(T, RT)

promote_rule_rational(::Type{PT}, ::Type{RationalPoly{S, T}}) where {PT<:APL, S, T} = RationalPoly{promote_rule(PT, S), promote_rule(T, termtype(T))}
promote_rule_rational(::Type{RationalPoly{S, T}}, ::Type{RationalPoly{U, V}}) where {S, T, U, V} = RationalPoly{promote_rule(S, U), promote_rule(T, V)}

Base.promote_rule(::Type{RS}, ::Type{RT}) where {RS<:RationalPoly, RT<:RationalPoly} = promote_rule_rational(RS, RT)
Base.promote_rule(::Type{PT}, ::Type{RT}) where {PT<:APL, RT<:RationalPoly} = promote_rule_rational(PT, RT)
Base.promote_rule(::Type{RT}, ::Type{PT}) where {PT<:APL, RT<:RationalPoly} = promote_rule_rational(PT, RT)

#promote_rule(::Type{Term{C, U}}, ::Type{RationalPoly{C, S, T}}) where {C, S, T, U} = RationalPoly{C, promote_type(U, S), T}
#promote_rule(::Type{RationalPoly{C, S, T}}, ::Type{Term{C, U}}) where {C, S, T, U} = RationalPoly{C, promote_type(U, S), T}
#promote_rule(::Type{Polynomial{C, U}}, ::Type{RationalPoly{C, S, T}}) where {C, S, T, U} = RationalPoly{C, promote_type(U, S), T}
#promote_rule(::Type{RationalPoly{C, S, T}}, ::Type{Polynomial{C, U}}) where {C, S, T, U} = RationalPoly{C, promote_type(U, S), T}

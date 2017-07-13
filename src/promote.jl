Base.promote_rule{PT<:APL, PS<:APL}(::Type{PT}, ::Type{PS}) = promote_rule(polynomialtype(PT), polynomialtype(PS))
Base.promote_rule{PT<:APL}(::Type{PT}, ::Type{PT}) = PT

promote_rule_constant{T, NT, DT}(::Type{T}, ::Type{RationalPoly{NT, DT}}) = RationalPoly{promote_rule(T, NT), promote_rule(DT, termtype(DT))}

Base.promote_rule{T, PT<:APL}(::Type{PT}, ::Type{T}) = promote_rule_constant(T, PT)
Base.promote_rule{T, PT<:APL}(::Type{T}, ::Type{PT}) = promote_rule_constant(T, PT)
Base.promote_rule{T, RT<:RationalPoly}(::Type{T}, ::Type{RT}) = promote_rule_constant(T, RT)
Base.promote_rule{T, RT<:RationalPoly}(::Type{RT}, ::Type{T}) = promote_rule_constant(T, RT)

promote_rule_rational{PT<:APL, S, T}(::Type{PT}, ::Type{RationalPoly{S, T}}) = RationalPoly{promote_rule(PT, S), promote_rule(T, termtype(T))}
promote_rule_rational{S, T, U, V}(::Type{RationalPoly{S, T}}, ::Type{RationalPoly{U, V}}) = RationalPoly{promote_rule(S, U), promote_rule(T, V)}

Base.promote_rule{RS<:RationalPoly, RT<:RationalPoly}(::Type{RS}, ::Type{RT}) = promote_rule_rational(RS, RT)
Base.promote_rule{PT<:APL, RT<:RationalPoly}(::Type{PT}, ::Type{RT}) = promote_rule_rational(PT, RT)
Base.promote_rule{PT<:APL, RT<:RationalPoly}(::Type{RT}, ::Type{PT}) = promote_rule_rational(PT, RT)

#promote_rule{C, S, T, U}(::Type{Term{C, U}}, ::Type{RationalPoly{C, S, T}}) = RationalPoly{C, promote_type(U, S), T}
#promote_rule{C, S, T, U}(::Type{RationalPoly{C, S, T}}, ::Type{Term{C, U}}) = RationalPoly{C, promote_type(U, S), T}
#promote_rule{C, S, T, U}(::Type{Polynomial{C, U}}, ::Type{RationalPoly{C, S, T}}) = RationalPoly{C, promote_type(U, S), T}
#promote_rule{C, S, T, U}(::Type{RationalPoly{C, S, T}}, ::Type{Polynomial{C, U}}) = RationalPoly{C, promote_type(U, S), T}

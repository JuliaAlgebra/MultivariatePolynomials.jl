import Base.promote_rule

function promote_rule{V, C, S, T}(::Type{V}, ::Type{RationalPoly{C, S, T}})
    U = promote_type(V, S)
    RationalPoly{C, U, T}
end
promote_rule{S,T<:RationalPoly}(::Type{T}, ::Type{S}) = promote_rule(S, T)
promote_rule{PT<:PolyType, C, S, T}(::Type{PT}, ::Type{RationalPoly{C, S, T}}) = RationalPoly{C, S, T}
promote_rule{PT<:PolyType,T<:RationalPoly}(::Type{T}, ::Type{PT}) = T
promote_rule{C, S, T, U}(::Type{Term{C, U}}, ::Type{RationalPoly{C, S, T}}) = RationalPoly{C, promote_type(U, S), T}
promote_rule{C, S, T, U}(::Type{RationalPoly{C, S, T}}, ::Type{Term{C, U}}) = RationalPoly{C, promote_type(U, S), T}
promote_rule{C, S, T, U}(::Type{Polynomial{C, U}}, ::Type{RationalPoly{C, S, T}}) = RationalPoly{C, promote_type(U, S), T}
promote_rule{C, S, T, U}(::Type{RationalPoly{C, S, T}}, ::Type{Polynomial{C, U}}) = RationalPoly{C, promote_type(U, S), T}
promote_rule{C, S, T, U, V}(::Type{RationalPoly{C, S, T}}, ::Type{RationalPoly{C, U, V}}) = RationalPoly{C, promote_type(S, U), promote_type(T, V)}

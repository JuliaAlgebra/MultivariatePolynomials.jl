import Base.promote_rule

# Promotion with PolyVar and Monomial
promote_rule{C}(::Type{Monomial{C}}, ::Type{PolyVar{C}}) = Monomial{C}
promote_rule{C}(::Type{PolyVar{C}}, ::Type{Monomial{C}}) = Monomial{C}
#promote_rule{S<:Union{Monomial, PolyVar}, T<:Union{Monomial, PolyVar}}(::Type{S}, ::Type{T}) = Monomial{iscomm{S}}

promote_rule{S,T<:TermType}(::Type{S}, ::Type{T}) = Polynomial{iscomm(T), promote_type(S, eltype(T))}
promote_rule{S,T<:TermType}(::Type{T}, ::Type{S}) = Polynomial{iscomm(T), promote_type(S, eltype(T))}
promote_rule{S,T<:PolyType}(::Type{S}, ::Type{T}) = Term{iscomm(T), promote_type(S, Int)}
promote_rule{S,T<:PolyType}(::Type{T}, ::Type{S}) = Term{iscomm(T), promote_type(S, Int)}

# Promotion with Term
promote_rule{C,S,T}(::Type{Term{C, S}}, ::Type{Term{C, T}}) = Term{C, promote_type(S, T)}
promote_rule{S,C,T}(::Type{S}, ::Type{Term{C, T}}) = Term{C, promote_type(S, T)}
promote_rule{S,C,T}(::Type{Term{C, T}}, ::Type{S}) = Term{C, promote_type(S, T)}
promote_rule{C,T}(::Type{Monomial{C}}, ::Type{Term{C, T}}) = Term{C, T}
promote_rule{C,T}(::Type{Term{C, T}}, ::Type{Monomial{C}}) = Term{C, T}
promote_rule{C,T}(::Type{PolyVar{C}}, ::Type{Term{C, T}}) = Term{C, T}
promote_rule{C,T}(::Type{Term{C, T}}, ::Type{PolyVar{C}}) = Term{C, T}

# Promotion with Polynomial
promote_rule{C, S, T}(::Type{Polynomial{C, S}}, ::Type{Term{C, T}}) = Polynomial{C, promote_type(S, T)}
promote_rule{C, S, T}(::Type{Term{C, T}}, ::Type{Polynomial{C, S}}) = Polynomial{C, promote_type(S, T)}

promote_rule{S<:TermType,T<:TermType}(::Type{S}, ::Type{T}) = Polynomial{iscomm(T), promote_type(eltype(S), eltype(T))}
promote_rule{S<:PolyType,T<:TermType}(::Type{T}, ::Type{S}) = Polynomial{iscomm(T), promote_type(Int, eltype(T))}
promote_rule{S<:PolyType,T<:TermType}(::Type{S}, ::Type{T}) = Polynomial{iscomm(T), promote_type(Int, eltype(T))}

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

# Hack see https://github.com/JuliaLang/julia/pull/18218
import Base.promote_op
promote_op{S<:PolyType,T}(*, ::Type{T}, ::Type{S}) = Any
promote_op{S<:PolyType,T}(*, ::Type{S}, ::Type{T}) = Any
promote_op{S<:PolyType,T<:PolyType}(*, ::Type{S}, ::Type{T}) = Any

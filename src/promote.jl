import Base.promote_rule
promote_rule{S,T<:TermType}(::Type{S}, ::Type{T}) = VecPolynomial{iscomm(T), promote_type(S, eltype(T))}
promote_rule{S,T<:TermType}(::Type{T}, ::Type{S}) = VecPolynomial{iscomm(T), promote_type(S, eltype(T))}
promote_rule{S,T<:PolyType}(::Type{S}, ::Type{T}) = Term{iscomm(T), promote_type(S, Int)}
promote_rule{S,T<:PolyType}(::Type{T}, ::Type{S}) = Term{iscomm(T), promote_type(S, Int)}

# Promotion with Term
promote_rule{S,C,T}(::Type{S}, ::Type{Term{C, T}}) = Term{C, promote_type(S, T)}
promote_rule{S,C,T}(::Type{Term{C, T}}, ::Type{S}) = Term{C, promote_type(S, T)}
promote_rule{C,T}(::Type{Monomial{C}}, ::Type{Term{C, T}}) = Term{C, T}
promote_rule{C,T}(::Type{Term{C, T}}, ::Type{Monomial{C}}) = Term{C, T}
promote_rule{C,T}(::Type{PolyVar{C}}, ::Type{Term{C, T}}) = Term{C, T}
promote_rule{C,T}(::Type{Term{C, T}}, ::Type{PolyVar{C}}) = Term{C, T}

promote_rule{S<:Union{Monomial,PolyVar},T<:Union{Monomial,PolyVar}}(::Type{S}, ::Type{T}) = Monomial{iscomm(S)}
promote_rule{S<:TermType,T<:TermType}(::Type{S}, ::Type{T}) = VecPolynomial{iscomm(T), promote_type(eltype(S), eltype(T))}
promote_rule{S<:PolyType,T<:TermType}(::Type{T}, ::Type{S}) = VecPolynomial{iscomm(T), promote_type(Int, eltype(T))}
promote_rule{S<:PolyType,T<:TermType}(::Type{S}, ::Type{T}) = VecPolynomial{iscomm(T), promote_type(Int, eltype(T))}

function promote_rule{S,T<:RationalPoly}(::Type{S}, ::Type{T})
    U = promote_type(S, T.parameters[2], T.parameters[3])
    RationalPoly{iscomm(T), U, U}
end
promote_rule{S,T<:RationalPoly}(::Type{T}, ::Type{S}) = promote_rule(S, T)
promote_rule{S<:PolyType,T<:RationalPoly}(::Type{S}, ::Type{T}) = T
promote_rule{S<:PolyType,T<:RationalPoly}(::Type{T}, ::Type{S}) = T
promote_rule{S<:TermType,T<:RationalPoly}(::Type{S}, ::Type{T}) = promote_rule(eltype(S), T)
promote_rule{S<:TermType,T<:RationalPoly}(::Type{T}, ::Type{S}) = promote_rule(eltype(S), T)

# Hack see https://github.com/JuliaLang/julia/pull/18218
import Base.promote_op
promote_op{S<:PolyType,T}(*, ::Type{T}, ::Type{S}) = Any
promote_op{S<:PolyType,T}(*, ::Type{S}, ::Type{T}) = Any
promote_op{S<:PolyType,T<:PolyType}(*, ::Type{S}, ::Type{T}) = Any

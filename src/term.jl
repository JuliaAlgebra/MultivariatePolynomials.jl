export Term

struct Term{CoeffType, M <: AbstractMonomial} <: AbstractTerm{CoeffType}
    coefficient::CoeffType
    monomial::M
end

coefficient(t::Term) = t.coefficient
monomial(t::Term) = t.monomial
termtype(::Type{<:Term{C,M}}, ::Type{T}) where {C,M,T} = Term{T,M}
monomialtype(::Type{<:Term{C,M}}) where {C,M} = M

(t::Term)(s...) = substitute(Eval(), t, s)

LinearAlgebra.adjoint(t::Term) = Term(adjoint(coefficient(t)), monomial(t))

Base.convert(::Type{Term{T, M}}, m::AbstractMonomialLike) where {T, M} = Term(one(T), convert(M, m))

Base.promote_rule(::Type{Term{C,M1} where {C}}, M2::Type{<:AbstractMonomialLike}) where {M1} = (Term{C,promote_type(M1, M2)} where {C})
Base.promote_rule(M1::Type{<:AbstractMonomialLike}, ::Type{Term{C,M2} where {C}}) where {M2} = (Term{C,promote_type(M1, M2)} where {C})
Base.promote_rule(::Type{Term{C,M1} where {C}}, ::Type{Term{T,M2}}) where {T,M1,M2} = (Term{C,promote_type(M1, M2)} where {C})
Base.promote_rule(::Type{Term{T,M2}}, ::Type{Term{C,M1} where {C}}) where {T,M1,M2} = (Term{C,promote_type(M1, M2)} where {C})

combine(t1::Term, t2::Term) = combine(promote(t1, t2)...)
combine(t1::T, t2::T) where {T <: Term} = Term(t1.coefficient + t2.coefficient, t1.monomial)
compare(t1::Term, t2::Term) = monomial(t1) > monomial(t2)
function MA.promote_operation(::typeof(combine), ::Type{Term{S, M1}}, ::Type{Term{T, M2}}) where {S, T, M1, M2}
    return Term{MA.promote_operation(+, S, T), promote_type(M1, M2)}
end

function MA.mutability(::Type{Term{C,T}}) where {C,T}
    if MA.mutability(C) isa MA.IsMutable && MA.mutability(T) isa MA.IsMutable
        return MA.IsMutable()
    else
        return MA.IsNotMutable()
    end
end

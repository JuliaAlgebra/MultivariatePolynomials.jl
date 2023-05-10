struct Term{CoeffType,M<:AbstractMonomial} <: AbstractTerm{CoeffType}
    coefficient::CoeffType
    monomial::M
end

coefficient(t::Term) = t.coefficient
monomial(t::Term) = t.monomial
term_type(::Type{<:Term{C,M}}, ::Type{T}) where {C,M,T} = Term{T,M}
monomial_type(::Type{<:Term{C,M}}) where {C,M} = M
function Base.copy(t::Term)
    return Term(copy(t.coefficient), copy(t.monomial))
end

(t::Term)(s...) = substitute(Eval(), t, s)

LinearAlgebra.adjoint(t::Term) = Term(adjoint(coefficient(t)), monomial(t))

function Base.convert(::Type{Term{T,M}}, m::AbstractMonomialLike) where {T,M}
    return Term(one(T), convert(M, m))
end
function convert_constant(::Type{Term{C,M} where C}, α) where {M}
    return convert(Term{typeof(α),M}, α)
end

function Base.promote_rule(
    ::Type{Term{C,M1} where {C}},
    M2::Type{<:AbstractMonomialLike},
) where {M1}
    return (Term{C,promote_type(M1, M2)} where {C})
end
function Base.promote_rule(
    M1::Type{<:AbstractMonomialLike},
    ::Type{Term{C,M2} where {C}},
) where {M2}
    return (Term{C,promote_type(M1, M2)} where {C})
end
function Base.promote_rule(
    ::Type{Term{C,M1} where {C}},
    ::Type{Term{T,M2}},
) where {T,M1,M2}
    return (Term{C,promote_type(M1, M2)} where {C})
end
function Base.promote_rule(
    ::Type{Term{T,M2}},
    ::Type{Term{C,M1} where {C}},
) where {T,M1,M2}
    return (Term{C,promote_type(M1, M2)} where {C})
end
promote_rule_constant(::Type{T}, TT::Type{Term{C,M} where C}) where {T,M} = Any

combine(t1::Term, t2::Term) = combine(promote(t1, t2)...)
function combine(t1::T, t2::T) where {T<:Term}
    return Term(t1.coefficient + t2.coefficient, t1.monomial)
end
compare(t1::Term, t2::Term) = monomial(t1) < monomial(t2)
function MA.promote_operation(
    ::typeof(combine),
    ::Type{Term{S,M1}},
    ::Type{Term{T,M2}},
) where {S,T,M1,M2}
    return Term{MA.promote_operation(+, S, T),promote_type(M1, M2)}
end

function MA.mutability(::Type{Term{C,M}}) where {C,M}
    if MA.mutability(C) isa MA.IsMutable && MA.mutability(M) isa MA.IsMutable
        return MA.IsMutable()
    else
        return MA.IsNotMutable()
    end
end

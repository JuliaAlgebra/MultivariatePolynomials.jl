"""
    struct Term{CoeffType,M<:AbstractMonomial} <: AbstractTerm{CoeffType}
        coefficient::CoeffType
        monomial::M
    end

A representation of the multiplication between a `coefficient` and a `monomial`.

!!! note
    The `coefficient` does not need to be a `Number`. It can be for instance a
    multivariate polynomial. When computing a multivariate `gcd`, it is
    actually reformulated as a univariate `gcd` in one of the variable with
    coefficients being multivariate polynomials in the other variables.
    To create such a term, use [`term`](@ref) instead of `*`.
    For instance, if `p` is a polynomial and `m` is a monomial, `p * m` will
    multiply each term of `p` with `m` but `term(p, m)` will create a term
    with `p` as coefficient and `m` as monomial.
"""
struct Term{CoeffType,M<:AbstractMonomial} <: AbstractTerm{CoeffType}
    coefficient::CoeffType
    monomial::M
end

coefficient(t::Term) = t.coefficient
monomial(t::Term) = t.monomial
term_type(::Type{<:Term{C,M}}, ::Type{T}) where {C,M,T} = Term{T,M}
monomial_type(::Type{<:Term{C,M}}) where {C,M} = M

(t::Term)(s...) = substitute(Eval(), t, s)

function Base.convert(::Type{Term{T,M}}, m::AbstractMonomialLike) where {T,M}
    return Term(one(T), convert(M, m))
end
function Base.convert(::Type{Term{T,M}}, t::AbstractTerm) where {T,M}
    return Term{T,M}(convert(T, coefficient(t)), convert(M, monomial(t)))
end
function Base.convert(::Type{Term{T,M}}, t::Term{T,M}) where {T,M}
    return t
end

function convert_constant(::Type{Term{C,M} where C}, α) where {M}
    return convert(Term{typeof(α),M}, α)
end

function map_coefficients!(f::F, t::Term; nonzero = false) where {F<:Function}
    t.coefficient = f(t.coefficient)
    return t
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

# `Base.power_by_squaring` calls `Base.copy` and we want
# `t^1` to be a mutable copy of `t` so `copy` needs to be
# the same as `mutable_copy`.
Base.copy(t::Term) = MA.mutable_copy(t)
function MA.mutable_copy(t::Term)
    return Term(
        MA.copy_if_mutable(coefficient(t)),
        MA.copy_if_mutable(monomial(t)),
    )
end

function MA.operate_to!(
    t::Term,
    ::typeof(*),
    t1::AbstractTermLike,
    t2::AbstractTermLike,
)
    MA.operate_to!(t.coefficient, *, coefficient(t1), coefficient(t2))
    MA.operate_to!(t.monomial, *, monomial(t1), monomial(t2))
    return t
end

function MA.operate!(::typeof(*), t1::Term, t2::AbstractTermLike)
    MA.operate!(*, t1.coefficient, coefficient(t2))
    MA.operate!(*, t1.monomial, monomial(t2))
    return t1
end

# Needed to resolve ambiguity with
# MA.operate!(::typeof(*), ::AbstractMonomial, ::AbstractMonomialLike)
function MA.operate!(::typeof(*), t1::Term, t2::AbstractMonomialLike)
    MA.operate!(*, t1.coefficient, coefficient(t2))
    MA.operate!(*, t1.monomial, monomial(t2))
    return t1
end

function MA.operate!(::typeof(one), t::Term)
    MA.operate!(one, t.coefficient)
    MA.operate!(constant_monomial, t.monomial)
    return t
end

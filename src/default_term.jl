"""
    Term{CoeffType,M}

Type alias for [`StarAlgebras.Term{CoeffType,M}`](@ref).
A representation of the multiplication between a `coefficient` and a monomial.

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
const Term{CoeffType,M} = SA.Term{CoeffType,M}

monomial_type(::Type{<:SA.Term{<:Any,M}}) where {M} = M

(t::SA.Term)(s...) = substitute(Eval(), t, s)

function Base.convert(::Type{SA.Term{T,M}}, m::AbstractMonomialLike) where {T,M}
    return SA.Term(one(T), convert(M, m))
end
function Base.convert(::Type{SA.Term{T,M}}, t::SA.Term) where {T,M}
    return SA.Term{T,M}(convert(T, coefficient(t)), convert(M, monomial(t)))
end
function Base.convert(::Type{SA.Term{T,M}}, t::SA.Term{T,M}) where {T,M}
    return t
end
# Copy constructor needed for array operations (e.g., det)
SA.Term{T,M}(t::SA.Term{T,M}) where {T,M} = t

function convert_constant(::Type{SA.Term{C,M} where C}, α) where {M}
    return convert(SA.Term{typeof(α),M}, α)
end

function Base.promote_rule(
    ::Type{SA.Term{C,M1} where {C}},
    M2::Type{<:AbstractMonomialLike},
) where {M1}
    return (SA.Term{C,promote_type(M1, M2)} where {C})
end
function Base.promote_rule(
    M1::Type{<:AbstractMonomialLike},
    ::Type{SA.Term{C,M2} where {C}},
) where {M2}
    return (SA.Term{C,promote_type(M1, M2)} where {C})
end
function Base.promote_rule(
    ::Type{SA.Term{C,M1} where {C}},
    ::Type{SA.Term{T,M2}},
) where {T,M1,M2}
    return (SA.Term{C,promote_type(M1, M2)} where {C})
end
function Base.promote_rule(
    ::Type{SA.Term{T,M2}},
    ::Type{SA.Term{C,M1} where {C}},
) where {T,M1,M2}
    return (SA.Term{C,promote_type(M1, M2)} where {C})
end
promote_rule_constant(::Type{T}, TT::Type{SA.Term{C,M} where C}) where {T,M} = Any

combine(t1::SA.Term, t2::SA.Term) = combine(promote(t1, t2)...)
function combine(t1::T, t2::T) where {T<:SA.Term}
    return SA.Term(coefficient(t1) + coefficient(t2), monomial(t1))
end
function MA.promote_operation(
    ::typeof(combine),
    ::Type{SA.Term{S,M1}},
    ::Type{SA.Term{T,M2}},
) where {S,T,M1,M2}
    return SA.Term{MA.promote_operation(+, S, T),promote_type(M1, M2)}
end

# MA.mutability, Base.copy, and MA.mutable_copy for SA.Term are defined in StarAlgebras

function MA.operate_to!(
    t::SA.Term,
    ::typeof(*),
    t1::AbstractTermLike,
    t2::AbstractTermLike,
)
    MA.operate_to!(t.coefficient, *, coefficient(t1), coefficient(t2))
    MA.operate_to!(t.basis_element, *, monomial(t1), monomial(t2))
    return t
end

function MA.operate!(::typeof(*), t1::SA.Term, t2::AbstractTermLike)
    MA.operate!(*, t1.coefficient, coefficient(t2))
    MA.operate!(*, t1.basis_element, monomial(t2))
    return t1
end

# Needed to resolve ambiguity with
# MA.operate!(::typeof(*), ::AbstractMonomial, ::AbstractMonomialLike)
function MA.operate!(::typeof(*), t1::SA.Term, t2::AbstractMonomialLike)
    MA.operate!(*, t1.coefficient, coefficient(t2))
    MA.operate!(*, t1.basis_element, monomial(t2))
    return t1
end

function MA.operate!(::typeof(one), t::SA.Term)
    MA.operate!(one, t.coefficient)
    MA.operate!(constant_monomial, t.basis_element)
    return t
end

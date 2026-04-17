# SA.Term now has 3 type params: Term{T,A,I}
# T = coefficient type, A = algebra type, I = index type.
# The monomial is reconstructed via basis(algebra)[index].

# Convenience constructor: create a Term from coefficient + monomial
function SA.Term(coeff, mono::AbstractMonomialLike)
    alg = algebra(FullBasis{Monomial}(mono))
    idx = exponents(mono)
    return SA.Term(alg, idx, coeff)
end

# monomial_type: derive from the algebra
monomial_type(::Type{<:SA.Term{<:Any,A}}) where {A} = monomial_type(A)

(t::SA.Term)(s...) = substitute(Eval(), t, s)

# combine: same monomial (same index), add coefficients
function combine(t1::T, t2::T) where {T<:SA.Term}
    return SA.Term(t1.algebra, t1.index, coefficient(t1) + coefficient(t2))
end
function combine(t1::SA.Term, t2::SA.Term)
    # Different types: promote via algebra_element and sum
    return algebra_element(t1) + algebra_element(t2)
end

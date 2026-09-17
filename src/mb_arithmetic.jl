# SA.Term constructor for Polynomial{...} basis elements (monomials)
function SA.Term(coeff, p::Polynomial{B}) where {B}
    alg = algebra(FullBasis{B}(p))
    idx = exponents(p)
    return SA.Term(alg, idx, coeff)
end

function Base.zero(
    ::Type{SA.AlgebraElement{T,A,C}},
) where {T,A<:_PolynomialAlgebra,C}
    M = monomial_type(A)
    return zero(algebra_element(term(zero(T), constant_monomial(M))))
end

# Arithmetic for AlgebraElement polynomials.
# Since polynomials ARE AlgebraElements now, SA handles most arithmetic.

const _AE = AbstractPolynomial

# polynomial() on an AlgebraElement is the identity
polynomial(a::AbstractPolynomial) = a
polynomial(a::AbstractTermLike) = algebra_element(a)
function polynomial(a::AbstractPolynomialLike, ::Type{T}) where {T}
    return SA.AlgebraElement{T}(polynomial(a))
end

# _APL (monomials, terms) + _AE: convert to AE first.
# Note: SA already defines +(::Term, ::AlgebraElement); we only handle
# AbstractMonomialLike (not Term) to avoid ambiguity.
for op in [:+, :-, :*]
    @eval begin
        Base.$op(p::AbstractMonomialLike, q::_AE) = $op(algebra_element(p), q)
        Base.$op(p::_AE, q::AbstractMonomialLike) = $op(p, algebra_element(q))
        Base.$op(p::Polynomial{<:AbstractMonomialIndexed}, q::_AE) =
            $op(algebra_element(SA.Term(1, p)), q)
        Base.$op(p::_AE, q::Polynomial{<:AbstractMonomialIndexed}) =
            $op(p, algebra_element(SA.Term(1, q)))
    end
end

# Scalars use the polynomial's coefficient type.
for op in [:+, :-]
    @eval begin
        function Base.$op(p::Union{T,Number}, q::_AE{T}) where {T}
            i = implicit(q)
            return $op(constant_algebra_element(SA.basis(i), convert(T, p)), i)
        end
        function Base.$op(p::_AE{T}, q::Union{T,Number}) where {T}
            i = implicit(p)
            return $op(i, constant_algebra_element(SA.basis(i), convert(T, q)))
        end
    end
end

function term_element(α, p::Polynomial)
    return SA.algebra_element(SA.Term(α, p))
end
Base.:*(α::Number, p::Polynomial) = term_element(α, p)

# SA.Term constructor for Polynomial{...} basis elements (monomials)
function SA.Term(coeff, p::Polynomial{B}) where {B}
    alg = algebra(FullBasis{B}(p))
    idx = SA.basis(alg)[p]
    return SA.Term(alg, idx, coeff)
end

SA._algebra_operand(m::AbstractMonomialLike) = term(m)

function Base.zero(
    ::Type{SA.AlgebraElement{T,A,C}},
) where {
    T,
    B,
    V,
    E,
    A<:SA.StarAlgebra{
        <:Union{Variables,_NoncommutativeBasis},
        Polynomial{B,V,E},
    },
    C,
}
    M = monomial_type(V)
    return zero(T, algebra(FullBasis{B}(constant_monomial(M))))
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

# Bare monomials enter SA's term arithmetic.
for op in [:+, :-, :*]
    @eval begin
        Base.$op(p::AbstractMonomialLike, q::_AE) = $op(term(p), q)
        Base.$op(p::_AE, q::AbstractMonomialLike) = $op(p, term(q))
        Base.$op(p::Polynomial{<:AbstractMonomialIndexed}, q::_AE) =
            $op(SA.Term(1, p), q)
        Base.$op(p::_AE, q::Polynomial{<:AbstractMonomialIndexed}) =
            $op(p, SA.Term(1, q))
    end
end

# Extend finite bases before embedding scalars as constant terms.
for op in [:+, :-]
    @eval begin
        function Base.$op(
            p::Union{T,Number},
            q::Union{AbstractPolynomial{T},AbstractTerm{T}},
        ) where {T}
            i = implicit(q)
            return $op(SA._constant_term(i, p), i)
        end
        function Base.$op(
            p::Union{AbstractPolynomial{T},AbstractTerm{T}},
            q::Union{T,Number},
        ) where {T}
            i = implicit(p)
            return $op(i, SA._constant_term(i, q))
        end
    end
end

function term_element(α, p::Polynomial)
    return SA.algebra_element(SA.Term(α, p))
end
Base.:*(α::Number, p::Polynomial) = term_element(α, p)

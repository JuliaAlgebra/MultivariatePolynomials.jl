# SA.Term constructor for Polynomial{...} basis elements (monomials)
function SA.Term(coeff, p::Polynomial{<:AbstractMonomialIndexed})
    alg = algebra(FullBasis{Monomial}(p))
    idx = exponents(p)
    return SA.Term(alg, idx, coeff)
end

# Arithmetic for AlgebraElement polynomials.
# Since polynomials ARE AlgebraElements now, SA handles most arithmetic.

const _AE{T} = SA.AlgebraElement{T,<:SA.StarAlgebra{<:Variables}}

# polynomial() on an AlgebraElement is the identity
polynomial(a::SA.AlgebraElement) = a

# _APL (monomials, terms) + _AE: convert to AE first
for op in [:+, :-, :*]
    @eval begin
        Base.$op(p::AbstractTermLike, q::_AE) = $op(algebra_element(p), q)
        Base.$op(p::_AE, q::AbstractTermLike) = $op(p, algebra_element(q))
        # Polynomial{Monomial,...} basis elements also participate
        Base.$op(p::Polynomial{<:AbstractMonomialIndexed}, q::_AE) =
            $op(algebra_element(SA.Term(1, p)), q)
        Base.$op(p::_AE, q::Polynomial{<:AbstractMonomialIndexed}) =
            $op(p, algebra_element(SA.Term(1, q)))
    end
end

# _AE + _AE: delegate to SA's generic path with promote_bases
for op in [:+, :-, :*]
    @eval begin
        function Base.$op(p::_AE, q::_AE)
            _p, _q = SA.promote_bases(p, q)
            return MA.operate_to!(SA._preallocate_output($op, _p, _q), $op, _p, _q)
        end
    end
end

# _AE + scalar and scalar + _AE
for op in [:+, :-]
    @eval begin
        function Base.$op(p, q::_AE)
            i = implicit(q)
            return $op(constant_algebra_element(SA.basis(i), p), i)
        end
        function Base.$op(p::_AE, q)
            i = implicit(p)
            return $op(i, constant_algebra_element(SA.basis(i), q))
        end
    end
end

# Number * basis_element → term algebra element
function term_element(α, p::Polynomial{B}) where {B}
    return algebra_element(
        SA.SparseCoefficients((p.exponents,), (α,)),
        FullBasis{B}(variables(p)),
    )
end

Base.:*(α::Number, p::Polynomial) = term_element(α, p)

# TODO: implement substitute/subs/differentiate directly on AlgebraElement

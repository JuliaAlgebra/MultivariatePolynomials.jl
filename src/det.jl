# Scalar determinant, used for recursive computation of the determinant
LinearAlgebra.det(p::AbstractPolynomialLike{<:Number}) = p

# Matrix determinant by cofactor expansion, adapted from
# `LinearAlgebraX.cofactor_det`.
function det_impl(A::AbstractMatrix{T}) where {T}
    get_elem = let A = A
        (i, j) -> LinearAlgebra.det(A[i, j])
    end
    r = first(size(A))
    if 1 < r
        total = LinearAlgebra.det(zero(T))
        for i in Base.OneTo(r)
            a = get_elem(i, 1)
            if !iszero(a)
                ii = Base.OneTo(r) .!= i
                jj = 2:r
                B = A[ii, jj]
                x = det_impl(B)
                x = MA.operate!!(*, x, a)
                if iseven(i)
                    x = MA.operate!!(-, x)
                end
                total = MA.operate!!(+, total, x)
            end
        end
        total
    elseif isone(r)
        MA.copy_if_mutable(get_elem(1, 1))
    else
        error("unexpected")
    end
end

collect_if_not_already_matrix(m::Matrix) = m
collect_if_not_already_matrix(m::AbstractMatrix) = collect(m)

function det_impl_outer(m::AbstractMatrix{T}) where {T}
    if 0 < LinearAlgebra.checksquare(m)
        det_impl(collect_if_not_already_matrix(m))
    else
        LinearAlgebra.det(one(T))
    end
end

# Determinants of narrow integer type: `LinearAlgebra` seems to
# promote these to `Float64` to prevent them from overflowing. We
# instead promote to `BigInt` to keep things exact. In the case of
# `Bool` we also need to promote for type stability.

const NarrowIntegerTypes =
    Union{Bool,UInt8,Int8,UInt16,Int16,UInt32,Int32,UInt64,Int64,UInt128,Int128}

const NarrowIntegerPolynomialLike =
    AbstractPolynomialLike{T} where {T<:NarrowIntegerTypes}

promote_if_narrow(m::AbstractMatrix{<:AbstractPolynomialLike}) = m

function promote_if_narrow(m::AbstractMatrix{<:NarrowIntegerPolynomialLike})
    return map((p -> polynomial(p, BigInt)), m)
end

# For type stability, we want to promote termlikes to polynomiallikes
# before attempting to calculate the determinant.
promote_if_termlike(m::AbstractMatrix{<:AbstractPolynomialLike}) = m
promote_if_termlike(m::AbstractMatrix{<:AbstractTermLike}) = map(polynomial, m)

promote_if_necessary(m) = promote_if_termlike(promote_if_narrow(m))

function LinearAlgebra.det(m::AbstractMatrix{<:AbstractPolynomialLike})
    return det_impl_outer(promote_if_necessary(m))
end

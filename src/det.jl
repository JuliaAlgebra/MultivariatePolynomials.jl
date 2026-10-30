# The determinant of terms and polynomials is implemented by StarAlgebras.
function LinearAlgebra.det(M::AbstractMatrix{<:AbstractMonomialLike})
    return LinearAlgebra.det(map(term, M))
end

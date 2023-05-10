function LinearAlgebra.det(M::Matrix{<:AbstractPolynomialLike})
    m = size(M)[1]
    if m > 2
        return sum(
            (-1)^(i - 1) * M[i, 1] * LinearAlgebra.det(M[1:end.!=i, 2:end]) for
            i in 1:m
        )
    elseif m == 2
        return M[1, 1] * M[2, 2] - M[2, 1] * M[1, 2]
    else
        return M[1, 1]
    end
end

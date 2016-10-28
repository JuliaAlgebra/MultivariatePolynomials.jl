using SumOfSquares

function coucou{T,S}(A::AbstractMatrix{T}, x::AbstractVector{S})
  @show T
  @show S
    TS = Base.promote_op(*,T,S)
    @show one(T)
    @show one(S)
    @show one(T) * one(S)
  @show TS
end

@polyvar x1 x2 x3 x4 x5
x = [x1, x2, x3, x4, x5]

# The matrix under consideration
J = [1 -1  1  1 -1;
    -1  1 -1  1  1;
     1 -1  1 -1  1;
     1  1 -1  1 -1;
    -1  1  1 -1  1]

xs = x.^2
xs = Array{Monomial}(xs)
@show typeof(J)
@show typeof(xs)
coucou(J, xs)

using MultivariatePolynomials
using DynamicPolynomials
#using TypedPolynomials
using BenchmarkTools

let
    n = 4
    using MultivariatePolynomials, DynamicPolynomials
    @polyvar x[1:3] 
    Ms = [randn(n, n) for i in 1:4]
    Ms = map(M -> M + transpose(M), Ms)
    M = sum(x[i] .* Ms[i] for i in 1:3) + Ms[4]
    det(M)
    @time det(M)
end

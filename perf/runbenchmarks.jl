using MultivariatePolynomials
using BenchmarkTools

#   let
#       using MultiPoly
#       x, y = generators(MPoly{Float64}, :x, :y)
#       p = x + y + 2x^2 + y^3
#
#       println(@benchmark(evaluate($p, 1.0, 2.0)))
#   end

let
    @polyvar x y
    p = x + y + 2x^2 + y^3

    vars = [x, y]
    vals = [1.0, 2.0]

    println(@benchmark(($p)($vals, $vars)))

#    vars_out_of_order = [vars[2], vars[1]]
#    println(@benchmark(($p)($vals, $vars_out_of_order)))
end

using MultivariatePolynomials
using DynamicPolynomials
#using TypedPolynomials
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

    varst = (x, y)
    varsv = [x, y]
    valst = (1.0, 2.0)
    valsv = [1.0, 2.0]

    display(@benchmark(($p)(variables($p) => $valst)))
    display(@benchmark(($p)(variables($p) => $valsv)))
    display(@benchmark(($p)($varst => $valst)))
    display(@benchmark(($p)($varsv => $valsv)))

#    vars_out_of_order = [vars[2], vars[1]]
#    println(@benchmark(($p)($vals, $vars_out_of_order)))
end

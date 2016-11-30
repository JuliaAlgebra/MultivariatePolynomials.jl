using MultivariatePolynomials
using BenchmarkTools

let
    @polyvar x y
    p = x + y + 2x^2 + y^3

    vars = [x, y]
    vals = [1.0, 2.0]

    println(@benchmark(($p)($vals, $vars)))

    vars_out_of_order = [vars[2], vars[1]]
    println(@benchmark(($p)($vals, $vars_out_of_order)))
end

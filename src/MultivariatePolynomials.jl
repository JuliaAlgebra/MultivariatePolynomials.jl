__precompile__()

module MultivariatePolynomials

import Base: *, +, -, /, ^, ==,
    promote_rule, convert, show, isless, size, getindex,
    one, zero, transpose, isapprox, @pure, dot, copy

include("types.jl")
include("operators.jl")
#include("show.jl")
#include("substitution.jl")

end # module

__precompile__()

module MultivariatePolynomials

import Base: show, length, getindex, vect, isless, isempty, start, done, next, convert, dot, copy, eltype, zero, one

include("mono.jl")
include("poly.jl")
include("rational.jl")
include("measure.jl")
include("exp.jl")
include("promote.jl")
include("comp.jl")
include("alg.jl")
include("diff.jl")
include("subs.jl")
include("algebraicset.jl")

include("show.jl")

end # module

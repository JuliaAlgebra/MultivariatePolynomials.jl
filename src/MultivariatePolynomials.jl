__precompile__()

module MultivariatePolynomials

import Base: show, length, getindex, vect, isless, isempty, start, done, next, convert, dot, copy, eltype, zero, one

abstract PolyType

include("cmono.jl")
include("ncmono.jl")
include("mono.jl")

include("poly.jl")
include("rational.jl")
include("measure.jl")
include("exp.jl")
include("promote.jl")

include("comp.jl")
include("ccomp.jl")

include("alg.jl")
include("calg.jl")
include("ncalg.jl")

include("diff.jl")
include("subs.jl")
include("algebraicset.jl")

include("show.jl")

end # module

#__precompile__()

module MultivariatePolynomials

import Base.show, Base.length, Base.getindex, Base.vect, Base.isless, Base.isempty, Base.start, Base.done, Base.next, Base.convert, Base.dot

include("mono.jl")
include("poly.jl")
include("rational.jl")
include("exp.jl")
include("promote.jl")
include("comp.jl")
include("alg.jl")
include("diff.jl")
include("subs.jl")
include("show.jl")

end # module

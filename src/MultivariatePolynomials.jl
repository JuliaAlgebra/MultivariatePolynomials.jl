__precompile__()

module MultivariatePolynomials

import Base: show, length, getindex, vect, isless, isempty, start, done, next, convert, dot, copy, eltype, zero, one

abstract PolyType{C}
export iscomm
iscomm{C}(::PolyType{C}) = C
zero{C}(::Type{PolyType{C}}) = zero(Polynomial{C, Int})
one{C}(::Type{PolyType{C}}) = one(Polynomial{C, Int})
zero{C}(p::PolyType{C}) = zero(PolyType{C})
one{C}(p::PolyType{C}) = one(PolyType{C})

include("mono.jl")
include("poly.jl")
include("rational.jl")
include("measure.jl")
include("exp.jl")
include("promote.jl")

include("comp.jl")

include("alg.jl")
include("calg.jl")
include("ncalg.jl")

include("diff.jl")
include("subs.jl")
include("algebraicset.jl")
include("norm.jl")

include("show.jl")

end # module

using MultivariatePolynomials
using FactCheck

#include("mono.jl")
include("poly.jl")
#include("rational.jl")
include("promote.jl")
include("comp.jl")
include("alg.jl")
include("diff.jl")
include("subs.jl")
#include("show.jl")

FactCheck.exitstatus()

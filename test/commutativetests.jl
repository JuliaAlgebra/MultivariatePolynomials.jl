include("mutable_arithmetics.jl")

include("zip.jl")
include("variable.jl")
include("monomial.jl")
include("term.jl")
include("monomial_vector.jl")
include("polynomial.jl")
include("det.jl")

include("rational.jl")
isdefined(Mod, Symbol("@complex_polyvar")) && include("complex.jl")

include("promote.jl")
include("hash.jl")
include("norm.jl")

include("algebra.jl")
include("comparison.jl")

include("substitution.jl")
include("differentiation.jl")
include("division.jl")
include("gcd.jl")
include("chain_rules.jl")

include("show.jl")

include("example1.jl")
include("example2.jl")

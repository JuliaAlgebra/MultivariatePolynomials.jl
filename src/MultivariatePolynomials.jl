module MultivariatePolynomials

import LinearAlgebra

import MutableArithmetics as MA

import StarAlgebras as SA

include("types.jl")

include("zip.jl")
include("lazy_iterators.jl")

include("variable.jl")
include("monomial.jl")
include("term.jl")
include("polynomial.jl")
include("monomial_vector.jl")

include("rational.jl")

include("show.jl")
include("hash.jl")

include("promote.jl")
include("conversion.jl")

include("complex.jl")
include("operators.jl")
include("comparison.jl")

include("substitution.jl")
include("differentiation.jl")
include("antidifferentiation.jl")
include("division.jl")
include("gcd.jl")
include("det.jl")


# Polynomial bases and algebra representation shared with MultivariateBases.
include("mb_interface.jl")
include("mb_variables.jl")
include("mb_polynomial.jl")
include("mb_bases.jl")
include("mb_mstructures.jl")
include("mb_monomial_basis.jl")
include("mb_algebra.jl")
include("mb_arithmetic.jl")

include("deprecate.jl")

const _EXCLUDE_SYMBOLS = [Symbol(@__MODULE__), :eval, :include]

for sym in names(@__MODULE__; all = true)
    sym_string = string(sym)
    if sym in _EXCLUDE_SYMBOLS ||
       startswith(sym_string, "_") ||
       startswith(sym_string, "@_")
        continue
    end
    if !(
        Base.isidentifier(sym) ||
        (startswith(sym_string, "@") && Base.isidentifier(sym_string[2:end]))
    )
        continue
    end
    @eval export $sym
end

end # module

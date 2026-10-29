using Test
using LinearAlgebra

using MultivariatePolynomials
const MP = MultivariatePolynomials

# The implementation of variables used by the tests can be selected with
# `Pkg.test(test_args = ["DynamicPolynomials"])`. It defaults to the one of
# `SimplePolynomials.jl` which does not need any other package.
const IMPLEMENTATION = isempty(ARGS) ? :SimplePolynomials : Symbol(only(ARGS))
if IMPLEMENTATION == :SimplePolynomials
    include("SimplePolynomials.jl")
else
    @eval import $IMPLEMENTATION
end
const Mod = getfield(@__MODULE__, IMPLEMENTATION)

include("utils.jl")

if IMPLEMENTATION == :SimplePolynomials
    # Tests that do not depend on the implementation
    include("comparison.jl")
    include("invalidations.jl")
end
include("scalar_arithmetic.jl")
include("extensions.jl")
if IMPLEMENTATION != :TypedPolynomials
    include("sa_polynomial.jl")
end
include("commutativetests.jl")
if isdefined(Mod, Symbol("@ncpolyvar"))
    include("noncommutativetests.jl")
end
if IMPLEMENTATION == :TypedPolynomials
    # Monomials of TypedPolynomials are `isbits` so they do not allocate
    include("allocations.jl")
end

if VERSION < v"1.14-"
    include("mooncake.jl")
end

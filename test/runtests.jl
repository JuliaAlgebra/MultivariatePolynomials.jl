using Test
using LinearAlgebra

using MultivariatePolynomials
const MP = MultivariatePolynomials

# Tests that do not need any specific polynomial implementation
include("comparison.jl")

include("utils.jl")

# Taken from JuMP/test/solvers.jl
function try_import(name::Symbol)
    try
        @eval import $name
        return true
    catch e
        return false
    end
end

if try_import(:DynamicPolynomials)
    Mod = DynamicPolynomials
    include("commutativetests.jl")
    include("noncommutativetests.jl")
end

if try_import(:TypedPolynomials)
    Mod = TypedPolynomials
    include("commutativetests.jl")
    include("allocations.jl")
end

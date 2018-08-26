using Test
using LinearAlgebra

using MultivariatePolynomials
const MP = MultivariatePolynomials

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

global implementation = :DynamicPolynomials
if try_import(:DynamicPolynomials)
    Mod = DynamicPolynomials
    include("commutativetests.jl")
    include("noncommutativetests.jl")
end

global implementation = :TypedPolynomials
if try_import(:TypedPolynomials)
    Mod = TypedPolynomials
    include("commutativetests.jl")
end

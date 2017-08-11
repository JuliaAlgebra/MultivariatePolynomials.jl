using Base.Test

using MultivariatePolynomials

import DynamicPolynomials
Mod = DynamicPolynomials
#include("commutativetests.jl")
#include("noncommutativetests.jl")

import TypedPolynomials
Mod = TypedPolynomials
#include("commutativetests.jl")

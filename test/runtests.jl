using Base.Test

using MultivariatePolynomials
const MP = MultivariatePolynomials

include("custom.jl")

import DynamicPolynomials
Mod = DynamicPolynomials
include("commutativetests.jl")
#include("noncommutativetests.jl")

import TypedPolynomials
Mod = TypedPolynomials
include("commutativetests.jl")

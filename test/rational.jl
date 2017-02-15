@testset "RationalPoly" begin
    @test RationalPoly{true, Int, Int}(1) == 1
    @test typeof(RationalPoly{false, Int, Int}(1)) == RationalPoly{false, Int, Int}
    @inferred RationalPoly{true, Int, Int}(1)
    @polyvar x
    @test typeof(1 / x) == RationalPoly{true, Int, Int}
    @test typeof(1.0 / x) == RationalPoly{true, Float64, Int}
    @inferred 1 / x
    @inferred 1.0 / x
    @test typeof([2.0x / x^2, (x+x) / (1 + 2x^2)]) == Vector{RationalPoly{true, Float64, Int}}
    @test 2 * (1/x * (1-x)) + (1/x * x) * (1/x^2 * x^2) - (1-x)/x == (1-x)/x + 1
    @test (1/x + 1/x) / 2 == ((1 / (x^2 - 1) + (x+1)) - (x+1)) * ((x^2 - 1) / x)
    @test typeof(zero(1/x)) == Term{true, Int}
    @test MultivariatePolynomials.iszero(zero(1/x))
    @test typeof(zero(RationalPoly{true, Float64, Int})) == Polynomial{true, Float64}
    @test MultivariatePolynomials.iszero(zero(RationalPoly{true, Float64, Int}))
    @test typeof(x / 2) == Term{true, Float64}
    @test typeof((x + x^2) / 3.0) == Polynomial{true, Float64}
end

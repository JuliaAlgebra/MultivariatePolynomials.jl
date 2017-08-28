const MP = MultivariatePolynomials
@testset "Monomial" begin
    Mod.@polyvar x

    @test zero(x^2) == 0
    @test zero(x^2) isa AbstractPolynomial{Int}
    @inferred zero(x^2)
    @test one(x^2) == 1
    @test one(x^2) isa AbstractMonomial
    @inferred one(x^2)
    Mod.@polyvar y[1:7]
    m = y[1] * y[3] * y[5] * y[7]
    @test issorted(variables(y[2] * m), rev=true)
    @test issorted(variables(m * y[4]), rev=true)
    @test issorted(variables(y[6] * m), rev=true)

    @test nvariables(monovec([x^2, prod(y[2:4])])) == 4

    @test nterms(x^2) == 1
    @test @inferred(terms(x^2)) == [x^2]

    @test degree(x * y[2]^2, x) == 1
    @test degree(x * y[2]^2, y[1]) == 0
    @test degree(x * y[2]^2, y[2]) == 2

    @test_throws ErrorException variable(x^2)
    @test_throws ErrorException variable(x*y[1])
    @test_throws ErrorException variable(constantmonomial(typeof(x)))

    @test variable(x^1) == x
    @test variable(x^1) isa AbstractVariable
    @test variable(x^2 + x - x^2) == x
    @test variable(x^2 + x - x^2) isa AbstractVariable
    @test variable(1.0x) == x
    @test variable(1.0x) isa AbstractVariable
    @test_throws ErrorException variable(x + 2x) == x

    @test monic(x^2) == x^2

    @test MP._div(2x^2*y[1]^3, x*y[1]^2) == 2x*y[1]
end

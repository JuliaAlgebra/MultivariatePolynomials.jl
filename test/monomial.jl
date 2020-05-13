const MP = MultivariatePolynomials
@testset "Monomial" begin
    Mod.@polyvar x

    @test zero(x^2) == 0
    @test zero(x^2) isa AbstractPolynomial{Int}
    @inferred zero(x^2)
    @test (@inferred one(x^2)) == 1
    @test one(x^2) isa AbstractMonomial
    @test (@inferred one(typeof(x^2))) == 1
    @test one(typeof(x^2)) isa AbstractMonomial

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

    @test_throws InexactError variable(x^2)
    @test_throws InexactError variable(x*y[1])
    @test_throws InexactError variable(constantmonomial(typeof(x)))

    @test x != constantmonomial(typeof(x))
    @test constantmonomial(typeof(x)) != x

    m = x^2
    @test x != m
    @test m != x
    typetests(m)
    typetests([x^2, x^3])
    @test (@inferred polynomial(m)) isa AbstractPolynomial{Int}
    @test (@inferred polynomial(m, Float64)) isa AbstractPolynomial{Float64}

    @test variable(x^1) == x
    @test variable(x^1) isa AbstractVariable
    @test variable(x^2 + x - x^2) == x
    @test variable(x^2 + x - x^2) isa AbstractVariable
    @test variable(1.0x) == x
    @test variable(1.0x) isa AbstractVariable
    @test_throws InexactError variable(x + 2x) == x

    @test monic(x^2) == x^2

    @test MP._div(2x^2*y[1]^3, x*y[1]^2) == 2x*y[1]

    @test transpose(x) == x
    @test adjoint(x) == x
    @test transpose(x^2) == x^2
    @test adjoint(x^2) == x^2

    @testset "Effective variables" begin
        T = variable_union_type(x)
        @test x isa T
        @test y[2] isa T
        @test T[x, y[2]] == @inferred effective_variables(x * y[2])
        @test T[x, y[2]] == @inferred effective_variables(y[2] * x)
        @test T[x] == @inferred effective_variables(x * y[2]^0)
        @test T[x] == @inferred effective_variables(y[2]^0 * x)
        @test T[y[2]] == @inferred effective_variables(x^0 * y[2])
        @test T[y[2]] == @inferred effective_variables(y[2] * x^0)
        @test T[x, y[2]] == @inferred effective_variables(y[3]^0 * x * y[2])
        @test T[x, y[2]] == @inferred effective_variables(y[2] * y[3]^0 * x)
        @test T[x, y[2]] == @inferred effective_variables(y[2] * x * y[3]^0)
    end
end

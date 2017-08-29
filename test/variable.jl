import MultivariatePolynomials: AbstractVariable, similarvariable, @similarvariable

@testset "Variable" begin
    @testset "polyvar macro index set" begin
        Mod.@polyvar x y z
        Mod.@polyvar x[1:3] y z[1:2]
        @test length(x) == 3
        @test length(z) == 2
        @test x[1] > x[2] > x[3] > y > z[1] > z[2]
    end
    @testset "PolyVar" begin
        Mod.@polyvar x
        @test 1 != x
        @test x != 0
        @test copy(x) == x
        @test nvariables(x) == 1
        @test !isapproxzero(x)
        @test !iszero(x)
        @test zero(x) == 0
        @test iszero(zero(x))
        @test zero(x) isa AbstractPolynomial{Int}
        @inferred zero(x)
        @test one(x) == 1
        @test one(x) isa AbstractMonomial
        @inferred one(x)

        typetests(x)
        @test (@inferred polynomial(x)) isa AbstractPolynomial{Int}
        @test (@inferred polynomial(x, Float64)) isa AbstractPolynomial{Float64}

        @test nterms(x) == 1
        @test @inferred(terms(x)) == [x]

        Mod.@polyvar y
        @test degree(x, x) == 1
        @test degree(x, y) == 0
        @test length(exponents(x)) == 1
        @test first(exponents(x)) == 1
        @test isconstant(x) == false

        @test divides(x, x) == true
        @test divides(x, y) == false
    end
    @testset "Create similar variable" begin
        Mod.@polyvar x y
        f = x^2 + y

        z = similarvariable(f, Val{:z})
        @test z isa AbstractVariable

        z = similarvariable(typeof(f), Val{:z})
        @test z isa AbstractVariable

        @inferred similarvariable(f, Val{:z})
        @inferred similarvariable(typeof(f), Val{:z})

        w = similarvariable(f, :w)
        @test w isa AbstractVariable

        @similarvariable f o
        @test o isa AbstractVariable

        m = @similarvariable f u
        @test m isa AbstractVariable
    end
end

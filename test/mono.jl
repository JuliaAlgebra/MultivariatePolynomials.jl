@testset "PolyVar and Monomial tests" begin
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
        @test zero(x) isa AbstractTerm{Int}
        @inferred zero(x)
        @test one(x) == 1
        @test one(x) isa AbstractMonomial
        @inferred one(x)

        Mod.@polyvar y
        @test MP.exponent(x, x) == 1
        @test MP.exponent(x, y) == 0
        @test length(MP.exponents(x)) == 1
        @test first(MP.exponents(x)) == 1
        @test isconstant(x) == false

        @test divides(x, x) == true
        @test divides(x, y) == false
    end
    @testset "Monomial" begin
        Mod.@polyvar x
        @test zero(x^2) == 0
        @test zero(x^2) isa AbstractTerm{Int}
        @inferred zero(x^2)
        @test one(x^2) == 1
        @test one(x^2) isa AbstractMonomial
        @inferred one(x^2)
        Mod.@polyvar y[1:7]
        m = y[1] * y[3] * y[5] * y[7]
        @test issorted(variables(y[2] * m), rev=true)
        @test issorted(variables(m * y[4]), rev=true)
        @test issorted(variables(y[6] * m), rev=true)

        @test MP.exponent(x * y[2]^2, x) == 1
        @test MP.exponent(x * y[2]^2, y[1]) == 0
        @test MP.exponent(x * y[2]^2, y[2]) == 2

        @test_throws ErrorException variable(x^2)
        @test_throws ErrorException variable(x*y[1])
        @test_throws ErrorException variable(constantmonomial(typeof(x)))
        @test variable(x^1) == x
        @test variable(x^1) isa AbstractVariable

        @test MP._div(2x^2*y[1]^3, x*y[1]^2) == 2x*y[1]
    end
    @testset "Monomial Vector" begin
        Mod.@polyvar x y
        @test x > y
        @test x^2 > y^2
        X = [x^2, x*y, y^2]
        for (i, m) in enumerate(monomials((x, y), 2))
            @test m == X[i]
        end
        @test length(monovec([y, x])) == 2
        X = monovec([x, 1, x*y])
        @test X[2:3][1] == x
        @test X[2:3][2] == 1
        @test monovec(X[[3, 2]])[1] == x
        @test monovec(X[[3, 2]])[2] == 1
        # Documentation examples
        @test monovec([x*y, x, x*y, x^2*y, x]) == [x^2*y, x*y, x]
        @test monovectype([x*y, x, 1, x^2*y, x]) <: AbstractVector{typeof(x*y)}
        @test monovectype([x*y, x, x*y, x^2*y, x]) <: AbstractVector
        σ, smv = sortmonovec([x*y, x, x*y, x^2*y, x])
        @test smv == [x^2*y, x*y, x]
        @test σ[1] == 4
        @test σ[2] in (1, 3)
        @test σ[3] in (2, 5)
        @test mergemonovec([[x*y, x, x*y], [x^2*y, x]]) == [x^2*y, x*y, x]
    end
end

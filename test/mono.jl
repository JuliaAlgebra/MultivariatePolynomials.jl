@testset "PolyVar and Monomial tests" begin
    @testset "polyvar macro index set" begin
        n = 3
        @polyvar x[1:n] y z[1:n-1]
        @test isa(x, Vector{PolyVar})
        @test isa(y, PolyVar)
        @test isa(z, Vector{PolyVar})
        @test length(x) == 3
        @test length(z) == 2
        @test x[1] > x[2] > x[3] > y > z[1] > z[2]
    end
    @testset "PolyVar" begin
        @polyvar x
        @test copy(x) == x
        @test nvars(x) == 1
        @test zero(x) == 0
        @test typeof(zero(x)) == VecPolynomial{Int}
        @inferred zero(x)
        @test one(x) == 1
        @test typeof(one(x)) == VecPolynomial{Int}
        @inferred one(x)
    end
    @testset "Monomial" begin
        @polyvar x
        @test_throws ArgumentError Monomial([x], [1,0])
        @test zero(x^2) == 0
        @test typeof(zero(x^2)) == VecPolynomial{Int}
        @inferred zero(x^2)
        @test one(x^2) == 1
        @test typeof(one(x^2)) == VecPolynomial{Int}
        @inferred one(x^2)
    end
    @testset "MonomialVector" begin
        @polyvar x y
        X = [x^2,x*y,y^2]
        for (i, m) in enumerate(monomials([x,y], 2))
            @test m == X[i]
        end
    end
end
module newmodule
    using Base.Test
    import MultivariatePolynomials
    @testset "Polyvar macro hygiene" begin
        # Verify that the @polyvar macro works when the package has been activated
        # with `import` instead of `using`.
        MultivariatePolynomials.@polyvar x y
        @test isa(x, MultivariatePolynomials.PolyVar)
        @test isa(y, MultivariatePolynomials.PolyVar)
    end
end

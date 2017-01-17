@testset "PolyVar and Monomial tests" begin
    @testset "polyvar macro index set" begin
        n = 3
        @polyvar x[1:n] y z[1:n-1]
        @test isa(x, Vector{PolyVar{true}})
        @test isa(y, PolyVar{true})
        @test isa(z, Vector{PolyVar{true}})
        @test length(x) == 3
        @test length(z) == 2
        @test x[1] > x[2] > x[3] > y > z[1] > z[2]
    end
    @testset "PolyVar" begin
        @test zero(PolyVar{true}) == 0
        @test one(PolyVar{false}) == 1
        @polyvar x
        @test copy(x) == x
        @test nvars(x) == 1
        @test zero(x) == 0
        @test typeof(zero(x)) == VecPolynomial{true, Int}
        @inferred zero(x)
        @test one(x) == 1
        @test typeof(one(x)) == VecPolynomial{true, Int}
        @inferred one(x)
    end
    @testset "Monomial" begin
        @test zero(PolyVar{false}) == 0
        @test one(PolyVar{true}) == 1
        @polyvar x
        @test_throws ArgumentError Monomial{true}([x], [1,0])
        @test zero(x^2) == 0
        @test typeof(zero(x^2)) == VecPolynomial{true, Int}
        @inferred zero(x^2)
        @test one(x^2) == 1
        @test typeof(one(x^2)) == VecPolynomial{true, Int}
        @inferred one(x^2)
        @polyvar y[1:7]
        m = y[1] * y[3] * y[5] * y[7]
        @test issorted(vars(y[2] * m), rev=true)
        @test issorted(vars(m * y[4]), rev=true)
        @test issorted(vars(y[6] * m), rev=true)
    end
    @testset "MonomialVector" begin
        @polyvar x y
        @test_throws ArgumentError MonomialVector{true}([x], [[1], [1,0]])
        X = [x^2,x*y,y^2]
        for (i, m) in enumerate(monomials([x,y], 2))
            @test m == X[i]
        end
        X = MonomialVector([x, 1, x*y])
        @test vars(X) == [x, y]
        @test X.Z == [[1, 1], [1, 0], [0, 0]]
        @test isa(MonomialVector{true}([1]), MonomialVector{true})
        @test isa(MonomialVector{false}([1]), MonomialVector{false})
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

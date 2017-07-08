@testset "PolyVar and Monomial tests" begin
    @testset "polyvar macro index set" begin
        n = 3
        @polyvar x[1:n] y z[1:n-1]
        @test length(x) == 3
        @test length(z) == 2
        @test x[1] > x[2] > x[3] > y > z[1] > z[2]
    end
    @testset "PolyVar" begin
        @polyvar x
        @test copy(x) == x
        @test nvars(x) == 1
        @test zero(x) == 0
        @inferred zero(x)
        @test one(x) == 1
        @inferred one(x)
    end
    @testset "Monomial" begin
        @polyvar x
        @test zero(x^2) == 0
        @inferred zero(x^2)
        @test one(x^2) == 1
        @inferred one(x^2)
        @polyvar y[1:7]
        m = y[1] * y[3] * y[5] * y[7]
        @test issorted(vars(y[2] * m), rev=true)
        @test issorted(vars(m * y[4]), rev=true)
        @test issorted(vars(y[6] * m), rev=true)
    end
    @testset "MonomialVector" begin
        @polyvar x y
        X = [x^2,x*y,y^2]
        for (i, m) in enumerate(monomials([x, y], 2))
            @test m == X[i]
        end
    end
end
module newmodule
    using Base.Test
    import DynamicPolynomials
    @testset "Polyvar macro hygiene" begin
        # Verify that the @polyvar macro works when the package has been activated
        # with `import` instead of `using`.
        DynamicPolynomials.@polyvar x y
        @test isa(x, DynamicPolynomials.PolyVar)
        @test isa(y, DynamicPolynomials.PolyVar)
    end
end

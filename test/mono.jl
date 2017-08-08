@testset "PolyVar and Monomial tests" begin
    @testset "polyvar macro index set" begin
        Mod.@polyvar x y z
        n = 3
        Mod.@polyvar x[1:n] y z[1:n-1]
        @test length(x) == 3
        @test length(z) == 2
        @test x[1] > x[2] > x[3] > y > z[1] > z[2]
    end
    @testset "PolyVar" begin
        Mod.@polyvar x
        @test copy(x) == x
        @test nvariables(x) == 1
        @test zero(x) == 0
        @test zero(x) isa AbstractTerm{Int}
        @inferred zero(x)
        @test one(x) == 1
        @test one(x) isa AbstractTerm{Int}
        @inferred one(x)
    end
    @testset "Monomial" begin
        Mod.@polyvar x
        @test zero(x^2) == 0
        @test zero(x^2) isa AbstractTerm{Int}
        @inferred zero(x^2)
        @test one(x^2) == 1
        @test one(x^2) isa AbstractTerm{Int}
        @inferred one(x^2)
        Mod.@polyvar y[1:7]
        m = y[1] * y[3] * y[5] * y[7]
        @test issorted(variables(y[2] * m), rev=true)
        @test issorted(variables(m * y[4]), rev=true)
        @test issorted(variables(y[6] * m), rev=true)
    end
    @testset "Monomial Vector" begin
        Mod.@polyvar x y
        X = [x^2, x*y, y^2]
        for (i, m) in enumerate(monomials([x, y], 2))
            @test m == X[i]
        end
        @test length(monovec([y, x])) == 2
        X = monovec([x, 1, x*y])
        @test X[2:3][1] == x
        @test X[2:3][2] == 1
        @test X[[3, 2]][1] == x
        @test X[[3, 2]][2] == 1
    end
end

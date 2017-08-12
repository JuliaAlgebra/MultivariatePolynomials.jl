@testset "Term and Polynomial tests" begin
    @testset "Term" begin
        Mod.@polyvar x
        @test Any(1x) == 1x
        @test one(1x) == one(1.0x) == 1
        @test zero(1x) == zero(1.0x) == 0
        @test nvariables(0.0x) == 1
        @test nvariables(1x) == 1
        #@inferred one(1x)
        @inferred zero(1x)
        #@inferred one(1.0x)
        @inferred zero(1.0x)

        #@test_throws InexactError Int(2x)
    end
    @testset "Polynomial" begin
        Mod.@polyvar x
        @test polynomial(1 + x) == 1 + x
        @test one(1 + x) == one(1.0 + x) == 1
        @test zero(1 + x) == zero(1.0 + x) == 0
        @test 1 != 1 + x
        @test 2x == x + x^2 + x - x^2
        @test x + x^2 - x^2 - x != 2x
        @test x^2 + x != x^2 + x + 1
        #@inferred one(1 + x)
        @inferred zero(1 + x)
        #@inferred one(1.0 + x)
        @inferred zero(1.0 + x)

        @test (1.0 + x) * x == x^2 + x
        @test constantterm(1, x) * (1 - x) == 1 - x
        @test promote_type(typeof(1-x), typeof(x)) <: AbstractPolynomial{Int}
        @test x != 1 - x

        Mod.@polyvar y

        @test maxdeg(x*y + 2 + x^2*y + x + y) == 3
        @test mindeg(x*y + 2 + x^2*y + x + y) == 0
        @test extdeg(x*y + 2 + x^2*y + x + y) == (0, 3)
        @test nvariables(x + y - x) == 2
        @test nvariables(x + x^2) == 1

        p = polynomial([4, 9], [x, x*x])
        @test coefficients(p) == [9, 4]
        @test monomials(p)[1] == x^2
        @test monomials(p)[2] == x
        @test p == dot([4, 9], [x, x*x])

        @inferred polynomial(i -> float(i), [x, x*x])
        @inferred polynomial(i -> 3 - float(i), monovec([x*x, x]))
        for p in (polynomial(i -> float(i), [x, x*x]),
                  polynomial(i -> 3 - float(i), monovec([x*x, x])))
            @test coefficients(p) == [2.0, 1.0]
            @test monomials(p) == monovec([x^2, x])
        end

        @test transpose(x + y) == x + y
    end
    @testset "Graded Lex Order" begin
        Mod.@polyvar x y z
        p = 3*y^2 + 2*y*x
        @test coefficients(p) == [2, 3]
        @test monomials(p) == monovec([x*y, y^2])
        # Examples from p. 59 of the 4th edition of "Ideals, Varieties, and Algorithms" of Cox, Little and O'Shea
        f = 4*x*y^2*z + 4*z^2 - 5*x^3 + 7*x^2*z^2
        @test coefficients(f) == [7, 4, -5, 4]
        @test monomials(f) == monovec([x^2*z^2, x*y^2*z, x^3, z^2])
    end
end

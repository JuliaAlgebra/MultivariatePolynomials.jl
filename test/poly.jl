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

        @test leadingterm(2x^2) == 2x^2
        @test nterms(2x^2) == 1

        @test term(x) isa AbstractTerm
        @test term(x^2) == x^2
        @test term(1x^2) isa AbstractTerm
        @test term(1x) == x

        Mod.@polyvar y
        @test MP.exponent(2x^2, x) == 2
        @test MP.exponent(2x^2, y) == 0
        @test MP.exponent(2x^2, y) == 0

        @test_throws InexactError push!([1], 2x)
        @test_throws ErrorException push!([x^2], 2x)
    end
    @testset "Polynomial" begin
        Mod.@polyvar x

        @test terms(polynomial([1, x^2, x, 2x^2])) == [3x^2, x, 1]
        @test terms(polynomial([x^3, 2x^3, x^2, -2x^2, x^2, x, 2, -2], MP.SortedState())) == [3x^3, x]

        @test polynomial(1 + x) == 1 + x
        @test leadingterm(1 + x) == x
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

        @test terms(polynomial(1 + x + x^2 - x + x^2)) == [2x^2, 1]

        @test (1.0 + x) * x == x^2 + x
        @test constantterm(1, x) * (1 - x) == 1 - x
        @test promote_type(typeof(1-x), typeof(x)) <: AbstractPolynomial{Int}
        @test x != 1 - x

        @test term(x + x^2 - x) isa AbstractTerm
        @test term(x + x^2 - x) == x^2
        @test term(x - x) isa AbstractTerm
        @test iszero(term(x - x))
        @test_throws ErrorException term(x + x^2)

        Mod.@polyvar y

        @test maxdeg(x*y + 2 + x^2*y + x + y) == 3
        @test mindeg(x*y + 2 + x^2*y + x + y) == 0
        @test extdeg(x*y + 2 + x^2*y + x + y) == (0, 3)
        @test leadingterm(x*y + 2 + x^2*y + x + y) == x^2*y
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

        @test removemonomials(4x^2*y + x*y + 2x, [x*y]) == 4x^2*y + 2x

        @test_throws InexactError push!([1], x+1)
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

@testset "Comparison" begin
    @testset "Graded Lex Order" begin
        Mod.@polyvar x y z
        @test x > y > z
        @test x^2 * y > y^3 > z
        @test y^2 >= x
        # Examples from p. 58, 59 of the 4th edition of "Ideals, Varieties, and Algorithms" of Cox, Little and O'Shea
        @test x^1 * y^2 * z^3 > x^3 * y^2
        @test !(x^1 * y^2 * z^3 < x^3 * y^2)
        @test x^1 * y^2 * z^4 > x^1 * y^1 * z^5
        @test !(x^1 * y^2 * z^4 < x^1 * y^1 * z^5)
        @test x^4 * y^7 * z > x^4 * y^2 * z^3
        @test !(x^4 * y^7 * z < x^4 * y^2 * z^3)
        @test x * y^5 * z^2 < x^4 * y * z^3
        @test !(x * y^5 * z^2 > x^4 * y * z^3)
        @test x^5 * y * z > x^4 * y * z^2
        @test !(x^5 * y * z < x^4 * y * z^2)
    end
    @testset "Equality" begin
        @testset "Monomial equality" begin
            Mod.@polyvar x y
            @test 1 == constant_monomial(x)
            @test 2 != constant_monomial(x)
            @test 2 != x
            @test 2 != x * y
            @test x * y != x
            @test 1x * y == x * y
            @test 2x * y != x * y
            @test x == monomial(x)
            @test x * y^0 == x
            @test x != x^0 * y
            @test x^2 != x * y
        end
        @testset "Polynomial equality" begin
            Mod.@polyvar x y
            @test isapproxzero(0.0x, ztol = 0.0)
            @test polynomial(CustomTerms(x + 1 - x)) isa AbstractPolynomial
            @test MultivariatePolynomials.right_constant_eq(
                polynomial(CustomTerms(x + 1 - x)),
                1,
            )
            @test MultivariatePolynomials.right_constant_eq(
                CustomTerms(x + 1 - x),
                1,
            )
            @test 2 * CustomPoly(x + 1 - x) == 2
            @test 2 != CustomTerms(x + 1 - x)
            @test 3x^2 == CustomTerms(x - x + x^2) * 3
            @test CustomPoly(-x + x^2) != x^2
            @test 2 * x * y + 3 * y^2 == 3 * y^2 + 2 * y * x
            @test 3 * x * y + 2 * y^2 != 3 * y^2 + 2 * y * x
            @test x + y != x * (1 + y)
            @test x * y == 3x + 2x * y - x - x * y - 2x
            @test isapproxzero((1 + 1e-8)x - x, ztol = 1e-7)
            @test !isapproxzero((1 + 1e-6)x - x, ztol = 1e-7)
            @test isapprox(
                (2 - 1e-3) * x * y + (3 + 1e-3) * y^2,
                3 * y^2 + 2 * y * x,
                rtol = 1e-2,
            )
            @test !isapprox(
                (2 - 1e-3) * x * y + (3 + 1e-1) * y^2,
                3 * y^2 + 2 * y * x,
                rtol = 1e-2,
            )
            @test isapprox(
                1e-3 * x * y + 3 * y^2 + x^2,
                x^2 + 3 * y^2,
                rtol = 1e-2,
                ztol = 1e-2,
            )
            @test isapprox(
                1e-3 * x * y + 3 * y^2 + x^2,
                x^2 + 3 * y^2,
                rtol = 1e-2,
                atol = 1e-2,
            )
            @test isapprox(
                3 * y^2 + x^2,
                x^2 + 1e-3 * x * y + 3 * y^2,
                rtol = 1e-2,
                ztol = 1e-2,
            )
            @test isapprox(
                3 * y^2 + x^2,
                x^2 + 1e-3 * x * y + 3 * y^2,
                rtol = 1e-2,
                atol = 1e-2,
            )
            @test !isapprox(
                3 * y^2 + x^2,
                x^2 + 1e-1 * x * y + 3 * y^2,
                rtol = 1e-2,
                ztol = 1e-2,
            )
            @test !isapprox(
                3.0 * y^2 + x + x^2,
                x + 3 * y^2,
                rtol = 1e-2,
                ztol = 1e-2,
            )
            @test isapprox(x + 1 - x, 1 + 1e-8, rtol = 1e-7)
            @test !isapprox(1 + 1e-8, x + 1 - x, rtol = 1e-9)
        end
        @testset "RationalPoly equality" begin
            Mod.@polyvar x y
            @test (x^2 - x - 6) / (x + 2) != x + 3
            @test x - 3 == (x^2 - x - 6) / (x + 2)
            @test (x^2 - x - 6) / (x - 3) == x + 2
            @test 3 != 4x / 2x
            @test 4x / 2x == 2
            @test 3 != 4x / 2x
            @test 1 - 1 / x == (x - 1) / x
            @test 1 - 1 / x != (1 - x) / x
            @test x + x / x == 1 + x^2 / x
            @test x - x / x == -(1 - x^2 / x)
            @test (1 + x) / x - 1 == 1 / x
            @test isapprox((1 + 1e-8)x, (x * y) / y, rtol = 1e-7)
            @test isapproxzero(((1 + 1e-8)x - x) / y, ztol = 1e-7)
            @test !isapproxzero(((1 + 1e-8)x - y) / y, ztol = 1e-9)
            @test isapprox(((1 + 1e-8)x * y) / y^2, x / y, rtol = 1e-7)
            @test isapprox((2x) / x, 2.001, rtol = 1e-2)
            @test !isapprox(2.001, (2x) / x, rtol = 1e-4)
        end
    end
    @testset "Equality between a Polynomial and a type not defining zero #22" begin
        Mod.@polyvar x
        # Polynomial of multiple terms
        p = x + x^2
        @test nothing !== p
        @test p !== nothing
        @test Dict{Int,Int}() != p
        @test p !== Dict{Int,Int}()
        # Polynomial of one term
        p = x + x^2 - x
        @test p !== nothing
        @test p !== Dict{Int,Int}()
        # Polynomial of no term
        # See https://github.com/blegat/MultivariatePolynomials.jl/issues/22
        p = x - x
        @test p !== nothing
        @test p !== Dict{Int,Int}()
    end
    @testset "compare" begin
        lex = LexOrder
        grlex = Graded{lex}
        rinvlex = Reverse{InverseLexOrder}
        grevlex = Graded{rinvlex}
        Mod.@polyvar x y z
        # [CLO13, p. 58]
        @test compare(x * y^2 * z^3, x^3 * y^2, lex) < 0
        @test compare(x * y^2 * z^3, x^3 * y^2, grlex) > 0
        @test compare(x * y^2 * z^3, x^3 * y^2, rinvlex) < 0
        @test compare(x * y^2 * z^3, x^3 * y^2, grevlex) > 0
        @test compare(x * y^2 * z^4, x * y * z^5, lex) > 0
        @test compare(x * y^2 * z^4, x * y * z^5, grlex) > 0
        @test compare(x * y^2 * z^4, x * y * z^5, rinvlex) > 0
        @test compare(x * y^2 * z^4, x * y * z^5, grevlex) > 0
        # [CLO13, p. 59]
        @test compare(x^5 * y * z, x^4 * y * z^2, lex) > 0
        @test compare(x^5 * y * z, x^4 * y * z^2, grlex) > 0
        @test compare(x^5 * y * z, x^4 * y * z^2, rinvlex) > 0
        @test compare(x^5 * y * z, x^4 * y * z^2, grevlex) > 0
        # [CLO13] Cox, D., Little, J., & OShea, D.
        # *Ideals, varieties, and algorithms: an introduction to computational algebraic geometry and commutative algebra*.
        # Springer Science & Business Media, **2013**.
    end
end

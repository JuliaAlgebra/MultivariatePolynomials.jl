if !@isdefined(Mod)
    include("../SimplePolynomials.jl")
    import .SimplePolynomials as Mod
    include("../utils.jl")
end

@testset "Comparison" begin
    @testset "Monomial hash consistency" begin
        Mod.@polyvar x y
        for (a, b) in (
            (one(x), 1),
            (monomial(x), x),
            (monomial([x, y], [2, 0]), x^2),
            (monomial([x, y], [0, 3]), y^3),
            (monomial([x, y], [1, 1]), y * x),
        )
            @test isequal(a, b)
            for seed in (UInt(0), UInt(1), typemax(UInt))
                @test (@inferred hash(a, seed)) == hash(b, seed)
            end
            @test Dict(a => 42)[b] == 42
        end
    end
    @testset "Graded Lex Order" begin
        Mod.@polyvar x y z
        O = Graded(LexOrder())
        for p in (x, x * y, 2x, x + y, variables(x * y))
            @test (@inferred ordering(p)) === O
            @test (@inferred ordering(typeof(p))) === O
        end
        @test sort([(1, 0), (0, 2), (0, 1)]; lt = ordering(x)) ==
              [(0, 1), (1, 0), (0, 2)]
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
            @test 1 ≈ constant_monomial(x)
            @test 2 != constant_monomial(x)
            @test !(2 ≈ constant_monomial(x))
            @test 2 != x
            @test !(2 ≈ x)
            @test 2 != x * y
            @test !(2 ≈ x * y)
            @test x * y != x
            @test !(x * y ≈ x)
            @test 1x * y == x * y
            @test 1x * y ≈ x * y
            @test 2x * y != x * y
            @test !(2x * y ≈ x * y)
            @test x == monomial(x)
            @test x ≈ monomial(x)
            @test x * y^0 == x
            @test x * y^0 ≈ x
            @test x != x^0 * y
            @test !(x ≈ x^0 * y)
            @test x^2 != x * y
            @test !(x^2 ≈ x * y)
        end
        @testset "Polynomial equality" begin
            Mod.@polyvar x y
            @test isapproxzero(0.0x, ztol = 0.0)
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
            @test !isequal((x^2 - x - 6) / (x + 2), x + 3)
            @test x - 3 == (x^2 - x - 6) / (x + 2)
            @test isequal(x - 3, (x^2 - x - 6) / (x + 2))
            @test (x^2 - x - 6) / (x - 3) == x + 2
            @test isequal((x^2 - x - 6) / (x - 3), x + 2)
            @test 3 != 4x / 2x
            @test !isequal(3, 4x / 2x)
            @test 4x / 2x == 2
            @test isequal(4x / 2x, 2)
            @test 3 != 4x / 2x
            @test !isequal(3, 4x / 2x)
            @test 1 - 1 / x == (x - 1) / x
            @test isequal(1 - 1 / x, (x - 1) / x)
            @test 1 - 1 / x != (1 - x) / x
            @test !isequal(1 - 1 / x, (1 - x) / x)
            @test x + x / x == 1 + x^2 / x
            @test isequal(x + x / x, 1 + x^2 / x)
            @test x - x / x == -(1 - x^2 / x)
            @test isequal(x - x / x, -(1 - x^2 / x))
            @test (1 + x) / x - 1 == 1 / x
            @test isequal((1 + x) / x - 1, 1 / x)
            @test isapprox((1 + 1e-8)x, (x * y) / y, rtol = 1e-7)
            @test isapproxzero(((1 + 1e-8)x - x) / y, ztol = 1e-7)
            @test !isapproxzero(((1 + 1e-8)x - y) / y, ztol = 1e-9)
            @test isapprox(((1 + 1e-8)x * y) / y^2, x / y, rtol = 1e-7)
            @test isapprox((2x) / x, 2.001, rtol = 1e-2)
            @test !isapprox(2.001, (2x) / x, rtol = 1e-4)
        end
    end
    @testset "Equality with unrelated values #22" begin
        Mod.@polyvar x
        for p in (x, x^2, 2x, x + x^2, polynomial(x), zero(x), x / (x + 1))
            for other in (nothing, Dict{Int,Int}(), Int[], [1, 2], :x, "x")
                @test p != other
                @test other != p
                @test !isequal(p, other)
                @test !isequal(other, p)
            end
        end
    end
    lex = LexOrder
    grlex = Graded{lex}
    rinvlex = Reverse{InverseLexOrder}
    grevlex = Graded{rinvlex}
    @static if Symbol(Mod) == :DynamicPolynomials
        @testset "cmp $M" for M in [lex, grlex, rinvlex, grevlex]
            Mod.@polyvar x y z monomial_order = M
            # [CLO13, p. 58]
            sgn = (M == lex || M == rinvlex) ? -1 : 1
            @test sgn * cmp(x * y^2 * z^3, x^3 * y^2) > 0
            @test cmp(x * y^2 * z^4, x * y * z^5) > 0
            # [CLO13, p. 59]
            @test cmp(x^5 * y * z, x^4 * y * z^2) > 0
            # [CLO13] Cox, D., Little, J., & OShea, D.
            # *Ideals, varieties, and algorithms: an introduction to computational algebraic geometry and commutative algebra*.
            # Springer Science & Business Media, **2013**.
        end
    end

    @testset "Comparison with NaN" begin
        Mod.@polyvar x y
        for p in (NaN * one(x), polynomial(NaN * one(x)), (NaN * x) / x)
            @test p != NaN
            @test NaN != p
            @test isequal(p, NaN)
            @test isequal(NaN, p)
        end
        @testset "$poly1 and $poly2" for (poly1, poly2) in [
            (NaN * x, NaN * x),
            (NaN * x + 0, NaN * x + 0),
            (NaN * x, NaN * x + 0),
            (NaN * x / y, NaN * x / y),
            (x / (NaN * y), x / (NaN * y)),
            ((NaN * x) / (NaN * y), (NaN * x) / (NaN * y)),
            ((NaN * x + 0) / (NaN * y + 0), (NaN * x + 0) / (NaN * y + 0)),
        ]
            @test poly1 != poly2
            @test poly1 != poly1
            @test isequal(poly1, poly2)
            @test isequal(poly1, poly1)
            # RationalPoly equality multiplies and thus allocates
            if !(poly1 isa RationalPoly)
                @test (@allocated poly1 != poly2) == 0
                @test (@allocated poly1 != poly1) == 0
                @test (@allocated isequal(poly1, poly2)) == 0
                @test (@allocated isequal(poly1, poly1)) == 0
            end
        end
    end
end

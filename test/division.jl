const MP = MultivariatePolynomials
@testset "Division" begin
    @testset "GCD and LCM" begin
        Mod.@polyvar x y z
        @test gcd(x*y, y) == y
        @test lcm(x*y, y) == x*y
        @test gcd(x*y^2*z, y*z^3) == y*z
        @test lcm(x*y^2*z, y*z^3) == x*y^2*z^3
        @test gcd(x^2*y^7*z^3, x^4*y^5*z^2) == x^2*y^5*z^2
        @test lcm(x^2*y^7*z^3, x^4*y^5*z^2) == x^4*y^7*z^3
    end
    # Taken from
    # Ideals, Varieties, and Algorithms
    # Cox, Little and O'Shea, Fourth edition
    # They have been adapted to the grlex ordering
    @testset "Divides" begin
        Mod.@polyvar x y z
        @test divides(leadingmonomial(x+y), x) # leadingmonomial(x+y) will be x^1*y^0 -> tricky test !
        @test !divides(leadingmonomial(x^2+y), x)
        @test divides(x*y, x^2*y)
        @test divides(x*y, x*y^2)
        @test divides(y*z, x*y*z)
        @test !divides(y*z, x*z)
        @test !divides(x^2*y, x*y^2)
    end
    @testset "Leading function" begin
        Mod.@polyvar x y z
        @test @inferred(monic(2x^2 + 4y + 2)) == x^2 + 2y + 1
        # See page 60
        p = 4x*y^2*z + 4z^2 + 7x^2*z^2 - 5x^3
        @test leadingcoefficient(p) == 7
        @test leadingmonomial(p) == x^2*z^2
        @test leadingterm(p) == 7x^2*z^2
    end
    @testset "Division examples" begin
        Mod.@polyvar x y
        @test (@inferred div(x*y^2 + 1, x*y + 1)) == y
        @test (@inferred rem(x*y^2 + 1, x*y + 1)) == -y + 1
        @test (@inferred div(x*y^2 + x, y)) == x*y
        @test (@inferred rem(x*y^2 + x, y)) == x
        @test (@inferred rem(x^4 + x^3 + (1+1e-10)*x^2 + 1, x^2 + x + 1)) == 1
        @test (@inferred rem(x^4 + x^3 + (1+1e-10)*x^2 + 1, x^2 + x + 1; ztol=1e-11)) ≈ -((1+1e-10)-1)x + 1
    end
    @testset "Division by multiple polynomials examples" begin
        function testdiv(p, ps)
            q, r = @inferred divrem(p, ps)
            @test p == dot(q, ps) + r
            q, r
        end
        Mod.@polyvar x y
        # Example 1
        q, r = testdiv(x*y^2 + 1, [x*y + 1, y + 1])
        @test q == [y, -1]
        @test r == 2
        # Example 2
        q, r = testdiv(x^2*y + x*y^2 + y^2, [x*y - 1, y^2 - 1])
        @test q == [x + y, 1]
        @test r == x + y + 1
        # Example 4
        q, r = testdiv(x^2*y + x*y^2 + y^2, [y^2 - 1, x*y - 1])
        @test q == [x + 1, x]
        @test r == 2x + 1
        # Example 5
        q, r = testdiv(x*y^2 - x, [x*y - 1, y^2 - 1])
        @test q == [y, 0]
        @test r == -x + y
        q, r = testdiv(x*y^2 - x, [y^2 - 1, x*y - 1])
        @test q == [x, 0]
        @test r == 0

        @testset "Issue DynamicPolynomials#62" begin
            p = x^2 + x + 1
            q = rem(p, [x^2-y])
            @test q == x + y + 1
        end

        @test (@inferred rem(x^2*y + (1+1e-10)*x*y + 1, [x^2 + x, y + 1])) == 1
        @test (@inferred rem(x^2*y + (1+1e-10)*x*y + 1, [x^2 + x, y + 1]; ztol=1e-11)) == -((1+1e-10)-1)x + 1
    end
    @testset "Univariate gcd" begin
        Mod.@polyvar x
        @test -(x + 1) == @inferred gcd(x^2 - 1, x^2 + 2x + 1)
        @test -(x + 1) == @inferred gcd(x^2 + 2x + 1, x^2 - 1)
        @test -(x + 1) == @inferred gcd(x^2 - 1, x + 1)
        @test -(x + 1) == @inferred gcd(x + 1, x^2 - 1)
        @test   x + 1  == @inferred gcd(x - x, x + 1)
        @test   x + 1  == @inferred gcd(x + 1, x - x)
        @test       0  == @inferred gcd(x - x, x^2 - x^2)
        @test       0  == @inferred gcd(x^2 - x^2, x - x)
    end
    @testset "Multivariate gcd $T" for T in [Int, Rational{BigInt}]
        function _mult_test(a::Number, b)
            return @test iszero(maxdegree(b))
        end
        function _mult_test(a, b)
            @test iszero(rem(a, b))
            @test iszero(rem(b, a))
        end
        function mult_test(expected, a, b)
            g = @inferred gcd(a, b)
            @test g isa polynomialtype(a)
            _mult_test(expected, g)
        end
        function sym_test(a, b, g)
            mult_test(g, a, b)
            mult_test(g, b, a)
        end
        function triple_test(a, b, c)
            sym_test(a * c, b * c, gcd(a, b) * c)
            sym_test(b * a, c * a, gcd(b, c) * a)
            sym_test(c * b, a * b, gcd(c, a) * b)
        end
        Mod.@polyvar x y z
        o = one(T)
        # Inspired from https://github.com/JuliaAlgebra/MultivariatePolynomials.jl/issues/160
        f1 = o * x * y + o * x
        f2 = o * y^2
        f3 = o * x
        sym_test(f1, f2, 1)
        sym_test(f2, f3, 1)
        sym_test(f3, f1, x)
        triple_test(f1, f2, f3)
        a = (o * x + o * y^2) * (o * z^3 + o * y^2 + o * x)
        b = (o * x + o * y + o * z) * (o * x^2 + o * y)
        c = (o * x + o * y + o * z) * (o * z^3 + o * y^2 + o * x)
        sym_test(a, b, 1)
        sym_test(b, c, x + y + z)
        sym_test(c, a, z^3 + y^2 + x)
        triple_test(a, b, c)
    end
    @testset "lcm" begin
        Mod.@polyvar x
        @test -(x-1) * (x+1)^2 == @inferred lcm(x^2 - 1, x^2 + 2x + 1)
        @test -(x-1) * (x+1)^2 == @inferred lcm(x^2 + 2x + 1, x^2 - 1)
        @test x^2 - 1 == @inferred lcm(x^2 - 1, x - 1)
        @test x^2 - 1 == @inferred lcm(x - 1, x^2 - 1)
        @test       0 == @inferred lcm(x - x, x + 1)
        @test       0 == @inferred lcm(x + 1, x - x)
        @test       0 == @inferred lcm(x - x, x^2 - x^2)
        @test       0 == @inferred lcm(x^2 - x^2, x - x)
    end
end

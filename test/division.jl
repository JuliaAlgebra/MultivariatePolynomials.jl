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
        @test iszero(@inferred MP.proddiff(2x*y, 3y^2*x))
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

        @test (@inferred rem(x^2*y + (1+1e-10)*x*y + 1, [x^2 + x, y + 1])) == 1
        @test (@inferred rem(x^2*y + (1+1e-10)*x*y + 1, [x^2 + x, y + 1]; ztol=1e-11)) == -((1+1e-10)-1)x + 1
    end
    @testset "Determinant" begin
        Mod.@polyvar x y
        @test 2x^3 - x^2 == @inferred det([2x^2 x; x x])
        @test x^3 - x^2 == @inferred det([x^2 x; x x])
        @test x^2 - y^2 == @inferred det([x y; y x])
        @test -2x^3 + 2x^2*y + 2x*y^2 - 4x^2 - 2x*y + 4x + 4y == @inferred det([x+1 x*y-x+2 0; 2x^2-3y^3 x-y 2; x+y x^2 0])
    end
end

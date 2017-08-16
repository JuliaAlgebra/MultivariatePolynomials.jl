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
        @test divides(x*y, x^2*y)
        @test divides(x*y, x*y^2)
        @test divides(y*z, x*y*z)
        @test !divides(y*z, x*z)
        @test !divides(x^2*y, x*y^2)
    end
    @testset "Leading function" begin
        Mod.@polyvar x y z
        # See page 60
        p = 4x*y^2*z + 4z^2 + 7x^2*z^2 - 5x^3
        @test leadingcoefficient(p) == 7
        @test leadingmonomial(p) == x^2*z^2
        @test leadingterm(p) == 7x^2*z^2
    end
    @testset "Division examples" begin
        Mod.@polyvar x y
        @test iszero(@inferred MP.proddiff(2x*y, 3y^2*x))
        @test div(x*y^2 + 1, x*y + 1) == y
        @test rem(x*y^2 + 1, x*y + 1) == -y + 1
        @test div(x*y^2 + x, y) == x*y
        @test rem(x*y^2 + x, y) == x
    end
end

@testset "Division" begin
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
        @test div(x*y^2 + 1, x*y + 1) == y
        @test rem(x*y^2 + 1, x*y + 1) == -y + 1
    end
end

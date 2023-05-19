import MultivariatePolynomials: pair_zip, tuple_zip

@testset "zip tests" begin
    @testset "pair_zip" begin
        Mod.@polyvar x y z

        @test @inferred(pair_zip((x, y), (1, 2))) == (x => 1, y => 2)
        @test_throws ArgumentError pair_zip((x, y, z), (1, 2))
        @test_throws ArgumentError pair_zip((x, y), (1, 2, 3))
        @test @inferred(pair_zip((1, :x, r"x"), ("w", 1.0, +))) ==
              (1 => "w", :x => 1.0, r"x" => +)
        @test @inferred(pair_zip((1, :x, r"x") => ("w", 1.0, +))) ==
              (1 => "w", :x => 1.0, r"x" => +)
    end

    @testset "tuple_zip" begin
        Mod.@polyvar x y z

        @test @inferred(tuple_zip((x, y), (1, 2))) == ((x, 1), (y, 2))
        @test_throws ArgumentError tuple_zip((x, y, z), (1, 2))
        @test_throws ArgumentError tuple_zip((x, y), (1, 2, 3))
        @test @inferred(tuple_zip((1, :x, r"x"), ("w", 1.0, +))) ==
              ((1, "w"), (:x, 1.0), (r"x", +))
    end
end

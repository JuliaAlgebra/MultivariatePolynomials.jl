import MultivariatePolynomials: pairzip, tuplezip

# Tests taken from TypedPolynomials
@testset "pairzip" begin
    @polyvar x y z

    @test @inferred(pairzip((x, y), (1, 2))) == (x=>1, y=>2)
    @test_throws ArgumentError pairzip((x, y, z), (1, 2))
    @test_throws ArgumentError pairzip((x, y), (1, 2, 3))
    @test @inferred(pairzip((1, :x, r"x"), ("w", 1.0, +))) == (1=>"w", :x=>1.0, r"x"=>+)
    @test @inferred(pairzip((1, :x, r"x")=>("w", 1.0, +))) == (1=>"w", :x=>1.0, r"x"=>+)
end

@testset "tuplezip" begin
    @polyvar x y z

    @test @inferred(tuplezip((x, y), (1, 2))) == ((x, 1), (y, 2))
    @test_throws ArgumentError tuplezip((x, y, z), (1, 2))
    @test_throws ArgumentError tuplezip((x, y), (1, 2, 3))
    @test @inferred(tuplezip((1, :x, r"x"), ("w", 1.0, +))) == ((1, "w"), (:x, 1.0), (r"x", +))
end

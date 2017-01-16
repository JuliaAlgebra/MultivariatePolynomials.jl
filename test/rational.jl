@testset "RationalPoly" begin
    @test RationalPoly{true, Int, Int}(1) == 1
    @test typeof(RationalPoly{false, Int, Int}(1)) == RationalPoly{false, Int, Int}
    @inferred RationalPoly{true, Int, Int}(1)
    @polyvar x
    @test typeof(1 / x) == RationalPoly{true, Int, Int}
    @inferred 1 / x
end

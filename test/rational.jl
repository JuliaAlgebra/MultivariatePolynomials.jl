@testset "RationalPoly" begin
    @test RationalPoly{Int,Int}(1) == 1
    @test typeof(RationalPoly{Int,Int}(1)) == RationalPoly{Int,Int}
    @inferred RationalPoly{Int,Int}(1)
    @polyvar x
    @test typeof(1 / x) == RationalPoly{Int, Int}
    @inferred 1 / x
end

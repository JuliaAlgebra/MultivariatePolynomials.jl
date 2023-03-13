@testset "RationalPoly" begin
    #@test RationalPoly{true, Int, Int}(1) == 1
    #@test typeof(RationalPoly{false, Int, Int}(1)) == RationalPoly{false, Int, Int}
    #@inferred RationalPoly{true, Int, Int}(1)
    Mod.@polyvar x
    @test 1 / x isa RationalPoly
    @test 1.0 / x isa RationalPoly
    @inferred 1 / x
    @inferred 1.0 / x
    @test eltype([2.0x / x^2, (x+x) / (1 + 2x^2)]) <: RationalPoly
    @test 2 * (1/x * (1-x)) + (1/x * x) * (1/x^2 * x^2) - (1-x)/x == (1-x)/x + 1
    @test (1/x + 1/x) / 2 == ((1 / (x^2 - 1) + (x+1)) - (x+1)) * ((x^2 - 1) / x)
    @test (x / (x + 1)) / (x - 1) == (x / (x - 1)) / (x + 1)
    #@test typeof(zero(1/x)) == Term{true, Int}
    @test iszero(zero(1/x))
    @test zero(1/x) == 0
    @test one(1/x) == 1
    #@test typeof(zero(RationalPoly{true, Float64, Int})) == Polynomial{true, Float64}
    #@test iszero(zero(RationalPoly{true, Float64, Int}))
    #@test typeof(x / 2) == Term{true, Float64}
    #@test typeof((x + x^2) / 3.0) == Polynomial{true, Float64}
    @test nothing !== x / (x^2 + 1)
    @test (x^2 + 1) / (2x) != nothing
    @test Dict{Int,Int}() != x / (x^2 + 1)
    @test (x^2 + 1) / (2x) != Dict{Int,Int}()
    @test numerator(x / x^2) == x
    @test denominator(x / x^2) == x^2
    @test inv(x / x^2) == x
    @test x / x^2 == inv(x)
    @test isone(((x+1) / (x-1)) / ((x+1) / (x-1)))
    @test ((x+1)^2 / (x-1)) / ((x+1) / (x-1)) == x+1
    poly = (x+1)/(x+2.0)
    RType = typeof(poly)
    @test RType(true) == one(RType)
    @test RType(false) == zero(RType)
    @test one(RType) isa RType
    @test zero(RType) isa RType
    @test one(poly) isa RType
    @test zero(poly) isa RType
end

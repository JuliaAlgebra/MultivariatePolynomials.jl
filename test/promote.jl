@testset "Promotion" begin
    @polyvar x y
    @inferred x*y+x
    @test typeof([x, x*y+x, x]) == Vector{VecPolynomial{true, Int}}
    @test typeof([1, x/y, x]) == Vector{RationalPoly{true, Int, Int}}
    @test typeof([(x^2-2x+2) x; x x^2]) == Matrix{VecPolynomial{true, Int}}
end

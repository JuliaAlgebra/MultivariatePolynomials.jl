@testset "Promotion" begin
    @polyvar x y
    @test typeof([x, x*y+x, x]) == Vector{VecPolynomial{Int}}
    @test typeof([1, x/y, x]) == Vector{RationalPoly{Int, Int}}
    @test typeof([(x^2-2x+2) x; x x^2]) == Matrix{VecPolynomial{Int}}
end

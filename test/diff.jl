@testset "Differentiation" begin
    @polyvar x y
    @test differentiate(MatPolynomial{Int}((i,j)->1, [x]), y) == 0
    @test differentiate(x*y + 3y^2 , [x, y]) == [y, x+6y]
end

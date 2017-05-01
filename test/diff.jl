@testset "Differentiation" begin
    @polyvar x y
    @test differentiate(MatPolynomial{true, Int}((i,j)->1, [x]), y) == 0
    @test differentiate(x*y + 3y^2 , [x, y]) == [y, x+6y]
    @test differentiate(1 / x , [x, y]) == [-1/x^2, 0]
    @test differentiate((x - y) / (x * y) , [x, y]) == [y^2 / (x * y)^2, -x^2 / (x * y)^2]
end

@testset "Differentiation" begin
    @polyvar x y
    @test differentiate(true*x+true*x^2, y) == 0
    @inferred differentiate(true*x+true*x^2, y)
    @test differentiate(MatPolynomial{true, Int}((i,j)->1, [x]), y) == 0
    @test differentiate(x*y + 3y^2 , [x, y]) == [y, x+6y]
    @inferred differentiate(x*y + 3y^2 , x)
    @inferred differentiate(x*y + 3y^2 , [x, y])
    @test differentiate(1 / x , [x, y]) == [-1/x^2, 0]
    @test differentiate((x - y) / (x * y) , [x, y]) == [y^2 / (x * y)^2, -x^2 / (x * y)^2]
    p = x^2-y+x*y
    @test differentiate(p, x, 2) == 2
    @test differentiate(p, x, 0) === p
    @test_throws DomainError differentiate(p, x, -1)
    @test differentiate(x, x, 1) == 1
    @inferred differentiate(x, x, 0)
    @inferred differentiate(x, x, 1)
    @test differentiate(x*y^4, y, 3) == 24x*y
    @inferred differentiate(x*y^4, y, 3)
    @test differentiate(x*y^4, y, 0) == x*y^4
    @test differentiate(2x^2, x, 2) == 4
    @inferred differentiate(2x^2, x, 2)
    @test differentiate(MatPolynomial{true, Int}((i,j)->1, [x, y]), y, 1) == 2x + 2y
    @inferred differentiate(MatPolynomial{true, Int}((i,j)->1, [x, y]), y, 0)
    @inferred differentiate(2x^2 + x*y + y^2, [x, y])
    p = differentiate(2x^2 + 3x*y + y^2, [x, y], 2)
    @test isa(p, Matrix{Polynomial{true,Int}})
    @test p == [4 3; 3 2]
    p = differentiate(2x^2 + 3x*y^2 + 4y^3 + 2.0, (x, y), 2)
    @test isa(p, Matrix{Polynomial{true,Float64}})
    @test p == [4.0 6.0y; 6.0y 6.0x+24.0y]
end

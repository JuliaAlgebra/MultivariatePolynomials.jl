@testset "Differentiation" begin
    @polyvar x y
    @test differentiate(true*x+true*x^2, y) == 0
    @inferred differentiate(true*x+true*x^2, y)
    @test differentiate(true*x+true*x^2, y) == 0
    @inferred differentiate(true*x+true*x^2, y)
    @test differentiate(MatPolynomial{true, Int}((i,j)->1, [x]), y) == 0
    @test differentiate(x*y + 3y^2 , [x, y]) == [y, x+6y]
    @test differentiate(x*y + 3y^2 , [x, y]) == [y, x+6y]
    @test differentiate(1 / x , [x, y]) == [-1/x^2, 0]
    @test differentiate((x - y) / (x * y) , [x, y]) == [y^2 / (x * y)^2, -x^2 / (x * y)^2]
    p = x^2-y+x*y
    @test differentiate(p, x, 2) == 2
    @test differentiate(p, x, 0) === p
    # -> Cannot `convert` an object of type String to an object of type DomainError
    # it seems DomainError does not accept any argument
    #@test_throws DomainError differentiate(p, x, -1)
    @test differentiate(x, x, 1) == 1
    @inferred differentiate(x, x, 0)
    @inferred differentiate(x, x, 1)
    @test differentiate(x*y^4, y, 3) == 24x*y
    @inferred differentiate(x*y^4, y, 3)
    @test differentiate(2x^2, x, 2) == 4
    @inferred differentiate(2x^2, x, 2)
    @test differentiate(MatPolynomial{true, Int}((i,j)->1, [x, y]), y, 1) == 2x + 2y
    @inferred differentiate(MatPolynomial{true, Int}((i,j)->1, [x, y]), y, 1)
end

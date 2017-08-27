@testset "Product with non-commutative monomials and variables" begin
    @ncpolyvar x y z
    a = x * (y * x)
    @test a.vars == [x, y, x]
    @test a.z == ones(Int, 3)
    a = y * (y * x)
    @test a.vars == [y, x]
    @test a.z == [2, 1]
    a = x*x*y^2*x
    @test a.vars == [x, y, x]
    @test a.z == [2, 2, 1]
    mv = MonomialVector([x, y, x*y])
    @test x * mv == MonomialVector([x^2, x*y, x^2*y])
    @test mv * x == MonomialVector([x^2, y*x, x*y*x])
    @test y^2 * mv == MonomialVector([y^2*x, y^3, y^2*x*y])
    @test mv * y^2 == MonomialVector([x*y^2, y^3, x*y^3])
    a = y*x^2
    @test a.vars == [y, x]
    @test a.z == [1, 2]
end

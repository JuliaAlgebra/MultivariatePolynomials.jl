@testset "Non-commutative MatPolynomial" begin
    @ncpolyvar x y
    P = MatPolynomial([2 3 4;
                       3 4 5;
                       4 5 6], [x*y, x^2, y^2])
    @test P.Q.Q == [4, 3, 5, 2, 4, 6]
    P = MatPolynomial((i,j) -> i + j, [x*y, x^2, y^2])
    @test P.Q.Q == [4, 3, 5, 2, 4, 6]
    p = polynomial(P)
    @test p.a == [4, 3, 5, 4, 3, 2, 6, 5, 4]
    @test p.x == monovec([x^4, x^3*y, x^2*y^2, x*y^3, x*y*x^2, x*y*x*y, y^4, y^2*x^2, y^2*x*y])
    @inferred MatPolynomial(Matrix{Float64}(0, 0), typeof(x)[]) == 0
    @test MatPolynomial(Matrix{Float64}(0, 0), typeof(x)[]) isa AbstractPolynomialLike{Float64}
    @test MatPolynomial(Matrix{Float64}(0, 0), typeof(x)[]) == 0
    P = MatPolynomial((i,j) -> ((i,j) == (1,1) ? 2 : 0), [x*y, x^2, y^2])
    Q = MatPolynomial([0 1; 1 0], [x^2, y^2])
    @test P != Q
end

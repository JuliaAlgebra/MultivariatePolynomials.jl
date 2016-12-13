@testset "Graded Lex Order" begin
    @polyvar x y z
    p = 3*y^2 + 2*y*x
    @test p.a == [2, 3]
    @test p.x == MonomialVector([x*y, y^2])
    # Examples from p. 59 of the 4th edition of "Ideals, Varieties, and Algorithms" of Cox, Little and O'Shea
    f = 4*x*y^2*z + 4*z^2 - 5*x^3 + 7*x^2*z^2
    @test f.a == [7, 4, -5, 4]
    @test f.x == MonomialVector([x^2*z^2, x*y^2*z, x^3, z^2])
end
@testset "MatPolynomial" begin
    @polyvar x y
    P = MatPolynomial{Int}((i,j) -> i + j, [x^2, x*y, y^2])
    p = VecPolynomial(P)
    @test p.a == [2, 6, 12, 10, 6]
    @test p.x == MonomialVector([x^4, x^3*y, x^2*y^2, x*y^3, y^4])
end

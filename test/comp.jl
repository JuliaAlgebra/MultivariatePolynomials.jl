@testset "Monomial equality" begin
    @polyvar x y
    @test x*y != x
    @test x == Monomial(x)
    @test Monomial([x,y], [1,0]) == x
    @test x != Monomial([x,y], [0,1])
    @test MonomialVector([x,y], [[1,0],[0,0]]) == MonomialVector([x], [[1],[0]])
end
@testset "Graded Lex Order" begin
    @polyvar x y z
    @test x > y > z
    @test x^2*y > y^3 > z
    # Examples from p. 58, 59 of the 4th edition of "Ideals, Varieties, and Algorithms" of Cox, Little and O'Shea
    @test x^1*y^2*z^3 > x^3*y^2
    @test !(x^1*y^2*z^3 < x^3*y^2)
    @test x^1*y^2*z^4 > x^1*y^1*z^5
    @test !(x^1*y^2*z^4 < x^1*y^1*z^5)
    @test x^4*y^7*z > x^4*y^2*z^3
    @test !(x^4*y^7*z < x^4*y^2*z^3)
    @test x*y^5*z^2 < x^4*y*z^3
    @test !(x*y^5*z^2 > x^4*y*z^3)
    @test x^5*y*z > x^4*y*z^2
    @test !(x^5*y*z < x^4*y*z^2)
end
@testset "Polynomial equality" begin
    @polyvar x y
    @test 2*x*y + 3*y^2 == 3*y^2 + 2*y*x
    @test 3*x*y + 2*y^2 != 3*y^2 + 2*y*x
    @test isapprox((2-1e-3)*x*y + (3+1e-3)*y^2, 3*y^2 + 2*y*x, rtol=1e-2)
    @test !isapprox((2-1e-3)*x*y + (3+1e-1)*y^2, 3*y^2 + 2*y*x, rtol=1e-2)
    @test isapprox(1e-3*x*y + 3*y^2 + x^2, x^2 + 3*y^2, rtol=1e-2, ztol=1e-2)
    @test isapprox(3*y^2 + x^2, x^2 + 1e-3*x*y + 3*y^2, rtol=1e-2, ztol=1e-2)
end

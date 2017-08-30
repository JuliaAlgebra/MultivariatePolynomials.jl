@testset "Product with non-commutative monomials and variables" begin
    Mod.@ncpolyvar x y z
    m1 = x * (y * x)
    @test variables(m1)[1] == x
    @test variables(m1)[2] == y
    @test variables(m1)[3] == x
    @test all(exponents(m1) .== 1)
    m2 = y * (y * x)
    @test variables(m2)[1] == y
    @test variables(m2)[2] == x
    @test collect(exponents(m2)) == [2, 1]
    for m3 in [(x*x)*(y^2*x), (x*x*y)*(y*x), x*(x*y^2*x), (x^2*y^2)*x]
        @test variables(m3)[1] == x
        @test variables(m3)[2] == y
        @test variables(m3)[3] == x
        @test collect(exponents(m3)) == [2, 2, 1]
    end
    mv = monovec([x, y, x*y])
    @test x * mv == [x^2, x*y, x^2*y]
    @test mv * x == [x^2, y*x, x*y*x]
    @test y^2 * mv == [y^2*x, y^3, y^2*x*y]
    @test mv * y^2 == [x*y^2, y^3, x*y^3]
    m4 = y*x^2
    @test variables(m4)[1] == y
    @test variables(m4)[2] == x
    @test collect(exponents(m4)) == [1, 2]

    @test (x*y) * (3x^2 + 2x*y + 4y^2) == 3x*y*x^2 + 2x*y*x*y + 4x*y^3

    P = [2 3 4;
         3 4 5;
         4 5 6]
    Q = [4 3 5;
         3 2 4
         5 4 6]
    p = @inferred polynomial(P, [x*y, x^2, y^2])
    @test coefficients(p) == [4, 3, 5, 4, 3, 2, 6, 5, 4]
    @test monomials(p) == [x^4, x^3*y, x^2*y^2, x*y^3, x*y*x^2, x*y*x*y, y^4, y^2*x^2, y^2*x*y]
    p = @inferred polynomial(Q, monovec([x*y, x^2, y^2]))
    @test coefficients(p) == [4, 3, 5, 4, 3, 2, 6, 5, 4]
    @test monomials(p) == [x^4, x^3*y, x^2*y^2, x*y^3, x*y*x^2, x*y*x*y, y^4, y^2*x^2, y^2*x*y]
    q = @inferred polynomial(Matrix{Float64}(0, 0), emptymonovec(x))
    @test q isa AbstractPolynomialLike{Float64}
    @test q == 0
end

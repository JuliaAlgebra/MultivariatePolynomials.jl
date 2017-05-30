@testset "Graded Lex Order" begin
    @polyvar x y z
    @test x > y > z
    @test x^2*y > y^3 > z
    @test y^2 >= x
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
@testset "Equality" begin
    @testset "Monomial equality" begin
        @polyvar x y
        @test 1 == Monomial([x, y], [0, 0])
        @test 2 != Monomial([x, y], [0, 0])
        @test 2 != x
        @test 2 != x*y
        @test 2 != MonomialVector([x, y], 1)
        @test x != MonomialVector([x, y], 1)
        @test x*y != x
        @test 1x*y == x*y
        @test 2x*y != x*y
        @test x == Monomial(x)
        @test Monomial([x, y], [1, 0]) == x
        @test x != Monomial([x, y], [0, 1])
        @test MonomialVector([x, y], [[1, 0], [0, 0]]) == MonomialVector([x], [[1], [0]])
        @test MonomialVector([x, y], 2) != MonomialVector([x, y], 1)
        z = x
        @polyvar x
        @test z != x
    end
    @testset "Polynomial equality" begin
        @polyvar x y
        @test 2*x*y + 3*y^2 == 3*y^2 + 2*y*x
        @test 3*x*y + 2*y^2 != 3*y^2 + 2*y*x
        @test x + y != x * (1 + y)
        @test x*y == 3x + 2x*y - x - x*y - 2x
        @test isapproxzero((1+1e-8)x - x, ztol=1e-7)
        @test !isapproxzero((1+1e-6)x - x, ztol=1e-7)
        @test isapprox((2-1e-3)*x*y + (3+1e-3)*y^2, 3*y^2 + 2*y*x, rtol=1e-2)
        @test !isapprox((2-1e-3)*x*y + (3+1e-1)*y^2, 3*y^2 + 2*y*x, rtol=1e-2)
        @test isapprox(1e-3*x*y + 3*y^2 + x^2, x^2 + 3*y^2, rtol=1e-2, ztol=1e-2)
        @test isapprox(3*y^2 + x^2, x^2 + 1e-3*x*y + 3*y^2, rtol=1e-2, ztol=1e-2)
        @test !isapprox(3*y^2 + x^2, x^2 + 1e-1*x*y + 3*y^2, rtol=1e-2, ztol=1e-2)
        @test !isapprox(3.0*y^2 + x + x^2, x + 3*y^2, rtol=1e-2, ztol=1e-2)
    end
    @testset "MatPolynomial equality" begin
        @polyvar x y
        @test MatPolynomial([2 3; 3 2], [x, y]) == 2x^2 + 2y^2 + 6x*y
    end
    @testset "SOSDecomposition equality" begin
        @polyvar x y
        @test !isapprox(SOSDecomposition{true}([x+y, x-y]), SOSDecomposition{true}([x+y]))
        @test !isapprox(SOSDecomposition{true}([x+y, x-y]), SOSDecomposition{true}([x+y, x+y]))
        @test isapprox(SOSDecomposition{true}([x+y, x-y]), SOSDecomposition{true}([x+y, x-y]))
        @test isapprox(SOSDecomposition{true}([x+y, x-y]), SOSDecomposition{true}([x-y, x+y+1e-8]), ztol=1e-7)
    end
    @testset "RationalPoly equality" begin
        @polyvar x y
        @test isapprox((1+1e-8)x, (x*y)/y, rtol=1e-7)
        @test isapproxzero(((1+1e-8)x - x)/y, ztol=1e-7)
        @test !isapproxzero(((1+1e-8)x - y)/y, ztol=1e-9)
        @test isapprox(((1+1e-8)x*y) / y^2, x / y, rtol=1e-7)
        @test isapprox((2x) / x, 2.001, rtol=1e-2)
        @test !isapprox(2.001,  (2x) / x, rtol=1e-4)
    end
    @testset "AtomicMeasure equality" begin
        @polyvar x y
        @test isapprox(AtomicMeasure([x, y], [1, 2], [[1, 2], [0, 1]]), AtomicMeasure([x, y], [2, 1], [[0, 1], [1, 2]]))
        @test !isapprox(AtomicMeasure([x, y], [1, 2], [[1, 2], [0, 1]]), AtomicMeasure([x, y], [1, 2], [[0, 1], [1, 2]]))
        @test !isapprox(AtomicMeasure([x, y], [1, 2], [[1, 2], [0, 1]]), AtomicMeasure([x, y], [1], [[1, 2]]))
    end
end
@testset "Equality between a Polynomial and a type not defining zero #22" begin
    @polyvar x
    # Polynomial of multiple terms
    p = x + x^2
    @test p != nothing
    @test p != Dict{Int,Int}()
    # Polynomial of one term
    p = x + x^2 - x
    @test p != nothing
    @test p != Dict{Int,Int}()
    # Polynomial of no term
    # See https://github.com/blegat/MultivariatePolynomials.jl/issues/22
    p = x - x
    @test p != nothing
    @test p != Dict{Int,Int}()
end

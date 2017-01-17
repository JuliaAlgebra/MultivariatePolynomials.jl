@testset "Term and Polynomial tests" begin
    @testset "TermContainer and TermType" begin
        @test eltype(MultivariatePolynomials.TermContainer{true, Int}) == Int
        @test eltype(MultivariatePolynomials.TermType{false, Float64}) == Float64
        @polyvar x
        @test eltype(MultivariatePolynomials.TermContainer{true}(x)) == Int
    end

    @testset "Term" begin
        @test eltype(Term{true, Int}) == Int
        @test zero(Term{false, Int}).α == 0
        @test one(Term{true, Int}).α == 1
        @polyvar x
        @test one(1x) == one(1.0x) == 1
        @test zero(1x) == zero(1.0x) == 0
        @test typeof(one(1x)) == Term{true, Int}
        @test typeof(zero(1x)) == Term{true, Int}
        @test typeof(one(1.0x)) == Term{true, Float64}
        @test typeof(zero(1.0x)) == Term{true, Float64}
        @inferred one(1x)
        @inferred zero(1x)
        @inferred one(1.0x)
        @inferred zero(1.0x)

        @test typeof(MultivariatePolynomials.TermContainer{true}(MultivariatePolynomials.TermContainer{true}(1))) == Term{true, Int}
        @inferred MultivariatePolynomials.TermContainer{true}(MultivariatePolynomials.TermContainer{true}(1))
    end
    @testset "VecPolynomial" begin
        @test eltype(VecPolynomial{true, Int}) == Int
        @polyvar x
        @test one(1 + x) == one(1.0 + x) == 1
        @test zero(1 + x) == zero(1.0 + x) == 0
        @test typeof(one(1 + x)) == VecPolynomial{true, Int}
        @test typeof(zero(1 + x)) == VecPolynomial{true, Int}
        @test typeof(one(1.0 + x)) == VecPolynomial{true, Float64}
        @test typeof(zero(1.0 + x)) == VecPolynomial{true, Float64}
        @inferred one(1 + x)
        @inferred zero(1 + x)
        @inferred one(1.0 + x)
        @inferred zero(1.0 + x)
        @polyvar y
        @test maxdeg(x*y + 2 + x^2*y + x + y) == 3
        @test mindeg(x*y + 2 + x^2*y + x + y) == 0
        @test extdeg(x*y + 2 + x^2*y + x + y) == (0, 3)
    end
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
        P = MatPolynomial{true, Int}((i,j) -> i + j, [x^2, x*y, y^2])
        p = VecPolynomial(P)
        @test p.a == [2, 6, 12, 10, 6]
        @test p.x == MonomialVector([x^4, x^3*y, x^2*y^2, x*y^3, y^4])
    end
    @testset "Non-commutative MatPolynomial" begin
        @ncpolyvar x y
        P = MatPolynomial{false, Int}((i,j) -> i + j, [x^2, x*y, y^2])
        p = VecPolynomial(P)
        @test p.a == [2, 3, 4, 5, 3, 4, 6, 4, 5]
        @test p.x == MonomialVector([x^4, x^3*y, x^2*y^2, x*y^3, x*y*x^2, x*y*x*y, y^4, y^2*x^2, y^2*x*y])
    end
end

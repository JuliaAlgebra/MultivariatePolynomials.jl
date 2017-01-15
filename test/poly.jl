@testset "Term and Polynomial tests" begin
    @testset "TermContainer and TermType" begin
        @test eltype(MultivariatePolynomials.TermContainer{Int}) == Int
        @test eltype(MultivariatePolynomials.TermType{Float64}) == Float64
    end

    @testset "Term" begin
        @test eltype(Term{Int}) == Int
        @test zero(Term{Int}).α == 0
        @test one(Term{Int}).α == 1
        @polyvar x
        @test one(1x) == one(1.0x) == 1
        @test zero(1x) == zero(1.0x) == 0
        @test typeof(one(1x)) == Term{Int}
        @test typeof(zero(1x)) == Term{Int}
        @test typeof(one(1.0x)) == Term{Float64}
        @test typeof(zero(1.0x)) == Term{Float64}
        @inferred one(1x)
        @inferred zero(1x)
        @inferred one(1.0x)
        @inferred zero(1.0x)

        @test typeof(MultivariatePolynomials.TermContainer(MultivariatePolynomials.TermContainer(1))) == Term{Int}
        @inferred MultivariatePolynomials.TermContainer(MultivariatePolynomials.TermContainer(1))
    end
    @testset "VecPolynomial" begin
        @test eltype(VecPolynomial{Int}) == Int
        @polyvar x
        @test one(1 + x) == one(1.0 + x) == 1
        @test zero(1 + x) == zero(1.0 + x) == 0
        @test typeof(one(1 + x)) == VecPolynomial{Int}
        @test typeof(zero(1 + x)) == VecPolynomial{Int}
        @test typeof(one(1.0 + x)) == VecPolynomial{Float64}
        @test typeof(zero(1.0 + x)) == VecPolynomial{Float64}
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
        P = MatPolynomial{Int}((i,j) -> i + j, [x^2, x*y, y^2])
        p = VecPolynomial(P)
        @test p.a == [2, 6, 12, 10, 6]
        @test p.x == MonomialVector([x^4, x^3*y, x^2*y^2, x*y^3, y^4])
    end
end

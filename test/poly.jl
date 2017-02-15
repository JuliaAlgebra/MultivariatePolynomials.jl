@testset "Term and Polynomial tests" begin
    @testset "TermContainer and TermType" begin
        @test zero(MultivariatePolynomials.TermContainer{false, Float64}) == 0
        @test eltype(MultivariatePolynomials.TermContainer{true, Int}) == Int
        @test eltype(MultivariatePolynomials.TermType{false, Float64}) == Float64
        @polyvar x
        @test eltype(MultivariatePolynomials.TermContainer(x)) == Int
    end

    @testset "Term" begin
        @test eltype(Term{true, Int}) == Int
        @test zero(Term{false, Int}).α == 0
        @test one(Term{true, Int}).α == 1
        @polyvar x
        @test typeof(Term(1x)) == Term{true, Int}
        @test Term(1x) == 1x
        @test typeof(Any(1x)) == Term{true, Int}
        @test Any(1x) == 1x
        @test one(1x) == one(1.0x) == 1
        @test zero(1x) == zero(1.0x) == 0
        @test typeof(one(1x)) == Term{true, Int}
        @test typeof(zero(1x)) == Term{true, Int}
        @test typeof(one(1.0x)) == Term{true, Float64}
        @test typeof(zero(1.0x)) == Term{true, Float64}
        @test eltype(1x) == Int
        @test eltype(1.0x^2) == Float64
        @test nvars(0.0x) == 1
        @test nvars(1x) == 1
        @inferred one(1x)
        @inferred zero(1x)
        @inferred one(1.0x)
        @inferred zero(1.0x)

        @test_throws InexactError Int(2x)

        @test typeof(MultivariatePolynomials.TermContainer(MultivariatePolynomials.TermContainer{true}(1))) == Term{true, Int}
        @inferred MultivariatePolynomials.TermContainer(MultivariatePolynomials.TermContainer{true}(1))
        @test !isempty(1x)
    end
    @testset "Polynomial" begin
        @test eltype(Polynomial{true, Int}) == Int
        @polyvar x
        @test_throws ArgumentError Polynomial{true, Int}([1, 2], [x])
        @test_throws ArgumentError Polynomial{true, Int}([1, 2], MonomialVector([x]))
        @test_throws InexactError Polynomial{true, Int}([1.5], [x])
        @test Polynomial(1 + x) == 1 + x
        @test one(1 + x) == one(1.0 + x) == 1
        @test zero(1 + x) == zero(1.0 + x) == 0
        @test typeof(one(1 + x)) == Polynomial{true, Int}
        @test typeof(zero(1 + x)) == Polynomial{true, Int}
        @test typeof(one(1.0 + x)) == Polynomial{true, Float64}
        @test typeof(zero(1.0 + x)) == Polynomial{true, Float64}
        @inferred one(1 + x)
        @inferred zero(1 + x)
        @inferred one(1.0 + x)
        @inferred zero(1.0 + x)
        @polyvar y
        @test maxdeg(x*y + 2 + x^2*y + x + y) == 3
        @test mindeg(x*y + 2 + x^2*y + x + y) == 0
        @test extdeg(x*y + 2 + x^2*y + x + y) == (0, 3)
        @test nvars(x + y - x) == 2
        @test nvars(x + x^2) == 1

        p = Polynomial([4, 9], [x, x*x])
        p.a == [9, 4]
        p.x[1] == x^2
        p.x[2] == x

        @inferred Polynomial(i -> float(i), [x, x*x])
        @inferred Polynomial(i -> 3 - float(i), MonomialVector([x*x, x]))
        for p in (Polynomial(i -> float(i), [x, x*x]),
                  Polynomial(i -> 3 - float(i), MonomialVector([x*x, x])))
            @test typeof(p) == Polynomial{true, Float64}
            @test p.a == [2.0, 1.0]
            @test p.x == MonomialVector([x^2, x])
        end

        @ncpolyvar ncpolyvar u v
        @inferred Polynomial(i -> i, [u, u*u, 1])
        p = Polynomial(i -> i, [u, u*u, 1])
        @test typeof(p) == Polynomial{false, Int}
        @test p.a == [2, 1, 3]
        @test p.x == MonomialVector([u^2, u, 1])

        @test u + v*u + 1 != v*u + u
        @test removemonomials(u + v*u + 1, [1, u*v]) == v*u + u
        @test removemonomials(u + u*v + 1, [u*v]) == 1 + u

        @inferred Polynomial(2u)
        @inferred Polynomial{false, Int}(2.0u)
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
        p = Polynomial(P)
        @test p.a == [2, 6, 12, 10, 6]
        @test p.x == MonomialVector([x^4, x^3*y, x^2*y^2, x*y^3, y^4])
        for i in 1:3
            for j in 1:3
                @test P[i, j] == i + j
            end
        end
        for P in (MatPolynomial((i,j) -> i * j, [y, x]),
                  MatPolynomial((i,j) -> (3-i) * (3-j), MonomialVector([y, x])),
                  MatPolynomial([1 2; 2 4], [y, x]),
                  MatPolynomial([4 2; 2 1], MonomialVector([y, x])))
            @test P.Q == [4, 2, 1]
            @test P.x[1] == x
            @test P.x[2] == y
        end
        P = MatPolynomial((i,j) -> ((i,j) == (1,1) ? 2 : 0), [x*y, x^2, y^2])
        Q = MatPolynomial([0 1; 1 0], [x^2, y^2])
        @test P == Q
    end
    @testset "Non-commutative MatPolynomial" begin
        @ncpolyvar x y
        P = MatPolynomial([2 3 4;
                           3 4 5;
                           4 5 6], [x*y, x^2, y^2])
        @test P.Q == [4, 3, 5, 2, 4, 6]
        P = MatPolynomial((i,j) -> i + j, [x*y, x^2, y^2])
        @test P.Q == [4, 3, 5, 2, 4, 6]
        p = Polynomial(P)
        @test p.a == [4, 3, 5, 4, 3, 2, 6, 5, 4]
        @test p.x == MonomialVector([x^4, x^3*y, x^2*y^2, x*y^3, x*y*x^2, x*y*x*y, y^4, y^2*x^2, y^2*x*y])
        @inferred MatPolynomial(Matrix{Float64}(), PolyVar{false}[]) == 0
        @test typeof(MatPolynomial(Matrix{Float64}(), PolyVar{false}[])) == MatPolynomial{false, Float64}
        @test MatPolynomial(Matrix{Float64}(), PolyVar{false}[]) == 0
        P = MatPolynomial((i,j) -> ((i,j) == (1,1) ? 2 : 0), [x*y, x^2, y^2])
        Q = MatPolynomial([0 1; 1 0], [x^2, y^2])
        @test P != Q
    end
    @testset "SOSDecomposition" begin
        @test isempty(SOSDecomposition(PolyVar{false}[]))
        @polyvar x y
        ps = [1, x + y, x^2, x*y, 1 + x + x^2]
        P = MatPolynomial(SOSDecomposition(ps))
        P.Q == [2 0 1 0 1; 0 1 0 0 0; 1 0 2 1 1; 0 0 1 1 0; 1 0 1 0 2]
        P.x == MonomialVector([x^2, x*y, x, y, 1])
        @test P == P
        @test isapprox(MatPolynomial(SOSDecomposition(P)), P)
    end
end

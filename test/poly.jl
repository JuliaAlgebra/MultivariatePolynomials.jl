@testset "Term and Polynomial tests" begin
    @testset "Term" begin
        @polyvar x
        @test Any(1x) == 1x
        @test one(1x) == one(1.0x) == 1
        @test zero(1x) == zero(1.0x) == 0
        @test nvariables(0.0x) == 1
        @test nvariables(1x) == 1
        #@inferred one(1x)
        @inferred zero(1x)
        #@inferred one(1.0x)
        @inferred zero(1.0x)

        #@test_throws InexactError Int(2x)
    end
    @testset "Polynomial" begin
        @polyvar x
        @test polynomial(1 + x) == 1 + x
        @test one(1 + x) == one(1.0 + x) == 1
        @test zero(1 + x) == zero(1.0 + x) == 0
        #@inferred one(1 + x)
        @inferred zero(1 + x)
        #@inferred one(1.0 + x)
        @inferred zero(1.0 + x)

        @test (1.0 + x) * x == x^2 + x
        @test constantterm(1, x) * (1 - x) == 1 - x
        @test promote_type(typeof(1-x), typeof(x)) <: AbstractPolynomial{Int}
        @test x != 1 - x

        @polyvar y

        @test maxdeg(x*y + 2 + x^2*y + x + y) == 3
        @test mindeg(x*y + 2 + x^2*y + x + y) == 0
        @test extdeg(x*y + 2 + x^2*y + x + y) == (0, 3)
        @test nvariables(x + y - x) == 2
        @test nvariables(x + x^2) == 1

        p = polynomial([4, 9], [x, x*x])
        p.a == [9, 4]
        p.x[1] == x^2
        p.x[2] == x
        @test p == dot([4, 9], [x, x*x])

        @inferred polynomial(i -> float(i), [x, x*x])
        @inferred polynomial(i -> 3 - float(i), monovec([x*x, x]))
        for p in (polynomial(i -> float(i), [x, x*x]),
                  polynomial(i -> 3 - float(i), monovec([x*x, x])))
            @test p.a == [2.0, 1.0]
            @test p.x == monovec([x^2, x])
        end

        @test transpose(x + y) == x + y
    end
    @testset "Graded Lex Order" begin
        @polyvar x y z
        p = 3*y^2 + 2*y*x
        @test coefficients(p) == [2, 3]
        @test monomials(p) == monovec([x*y, y^2])
        # Examples from p. 59 of the 4th edition of "Ideals, Varieties, and Algorithms" of Cox, Little and O'Shea
        f = 4*x*y^2*z + 4*z^2 - 5*x^3 + 7*x^2*z^2
        @test coefficients(f) == [7, 4, -5, 4]
        @test monomials(f) == monovec([x^2*z^2, x*y^2*z, x^3, z^2])
    end
    @testset "MatPolynomial" begin
        @polyvar x y
        P = MatPolynomial{Int}((i,j) -> i + j, [x^2, x*y, y^2])
        zP = zero(typeof(P))
        @test isempty(zP.Q)
        @test zP == 0
        p = polynomial(P)
        @test p.a == [2, 6, 12, 10, 6]
        @test p.x == monovec([x^4, x^3*y, x^2*y^2, x*y^3, y^4])
        for i in 1:3
            for j in 1:3
                @test P[i, j] == i + j
            end
        end
        for P in (MatPolynomial((i,j) -> i * j, [y, x]),
                  MatPolynomial((i,j) -> (3-i) * (3-j), monovec([y, x])),
                  MatPolynomial([1 2; 2 4], [y, x]),
                  MatPolynomial([4 2; 2 1], monovec([y, x])))
            @test P.Q.Q == [4, 2, 1]
            @test P.x[1] == x
            @test P.x[2] == y
        end
        P = MatPolynomial((i,j) -> ((i,j) == (1,1) ? 2 : 0), [x*y, x^2, y^2])
        Q = MatPolynomial([0 1; 1 0], [x^2, y^2])
        @test P == Q
        p = MatPolynomial([2 3; 3 2], [x, y])
        @test polynomial(p) isa AbstractPolynomial
        @test polynomial(p, Int) isa AbstractPolynomial
    end
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
#   @testset "SOSDecomposition" begin
#       @polyvar x y
#       @test isempty(SOSDecomposition(typeof(x)[]))
#       ps = [1, x + y, x^2, x*y, 1 + x + x^2]
#       P = MatPolynomial(SOSDecomposition(ps))
#       P.Q == [2 0 1 0 1; 0 1 0 0 0; 1 0 2 1 1; 0 0 1 1 0; 1 0 1 0 2]
#       P.x == [x^2, x*y, x, y, 1]
#       @test P == P
#       @test isapprox(MatPolynomial(SOSDecomposition(P)), P)
#   end
end

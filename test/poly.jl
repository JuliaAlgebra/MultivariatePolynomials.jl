@testset "Term and Polynomial tests" begin
    @testset "Term" begin
        @polyvar x
        @test Any(1x) == 1x
        @test_broken one(1x) == one(1.0x) == 1
        @test zero(1x) == zero(1.0x) == 0
        @test_broken nvars(0.0x) == 1
        @test_broken nvars(1x) == 1
        #@inferred one(1x)
        @inferred zero(1x)
        #@inferred one(1.0x)
        @inferred zero(1.0x)

        #@test_throws InexactError Int(2x)
    end
    @testset "Polynomial" begin
        @polyvar x
        @test_broken one(1 + x) == one(1.0 + x) == 1
        @test zero(1 + x) == zero(1.0 + x) == 0
        #@inferred one(1 + x)
        @inferred zero(1 + x)
        #@inferred one(1.0 + x)
        @inferred zero(1.0 + x)
        @polyvar y

        @test_broken maxdeg(x*y + 2 + x^2*y + x + y) == 3
        @test_broken mindeg(x*y + 2 + x^2*y + x + y) == 0
        @test_broken extdeg(x*y + 2 + x^2*y + x + y) == (0, 3)
        @test_broken nvars(x + y - x) == 2
        @test_broken nvars(x + x^2) == 1

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
    end
#   @testset "MatPolynomial" begin
#       @polyvar x y
#       P = MatPolynomial{true, Int}((i,j) -> i + j, [x^2, x*y, y^2])
#       p = Polynomial(P)
#       @test p.a == [2, 6, 12, 10, 6]
#       @test p.x == MonomialVector([x^4, x^3*y, x^2*y^2, x*y^3, y^4])
#       for i in 1:3
#           for j in 1:3
#               @test P[i, j] == i + j
#           end
#       end
#       for P in (MatPolynomial((i,j) -> i * j, [y, x]),
#                 MatPolynomial((i,j) -> (3-i) * (3-j), MonomialVector([y, x])),
#                 MatPolynomial([1 2; 2 4], [y, x]),
#                 MatPolynomial([4 2; 2 1], MonomialVector([y, x])))
#           @test P.Q == [4, 2, 1]
#           @test P.x[1] == x
#           @test P.x[2] == y
#       end
#       P = MatPolynomial((i,j) -> ((i,j) == (1,1) ? 2 : 0), [x*y, x^2, y^2])
#       Q = MatPolynomial([0 1; 1 0], [x^2, y^2])
#       @test P == Q
#   end
#   @testset "SOSDecomposition" begin
#       @polyvar x y
#       ps = [1, x + y, x^2, x*y, 1 + x + x^2]
#       P = MatPolynomial(SOSDecomposition(ps))
#       P.Q == [2 0 1 0 1; 0 1 0 0 0; 1 0 2 1 1; 0 0 1 1 0; 1 0 1 0 2]
#       P.x == [x^2, x*y, x, y, 1]
#       @test P == P
#       @test isapprox(MatPolynomial(SOSDecomposition(P)), P)
#   end
end

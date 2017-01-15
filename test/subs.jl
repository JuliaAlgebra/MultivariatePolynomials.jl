import Base.Test: @inferred

@testset "Substitution" begin
    @polyvar x[1:3]

    a = (x[1])([x[2]], [x[1]])
    b = x[2]
    @test (x[1])([x[2]], [x[1]]) == x[2]

    p = x[1]*x[2]*x[3]
    @test Int(p([1, 2, 3], x)) == 6

    p = x[1]^2 + x[1]*x[3] - 3
    @test Int(p([5, x[1]*x[2], 4], x)) == 42

    p = x[1]^2 + x[2]^2
    q = p([1 -1; 1 1] * x[1:2], x[1:2])
    @test q == 2p

    q = (x[1] + 1) / (x[1] + 2)
    @test isapproxzero(q([-1], [x[1]]))
    @test !isapproxzero(q([1], [x[1]]))
    @test isapprox(q([1], [x[1]]), 2/3)

    t = 2.0 * x[1] * x[2]^2
    @test t([1, 2, 3], x) == 8.0
    @test t([1, 3], [x[2], x[1]]) == 6.0

    P = [1 2 3; 2 4 5; 3 5 6]
    p = MatPolynomial(P, x)
    @test p(ones(3), x) == 31
    @test subs(p, ones(3), x) == 31

    p = x[1] + x[2] + 2*x[1]^2 + 3*x[1]*x[2]^2
    @inferred p([1.0, 2.0], [x[1], x[2]])

    @test subs(p, [2.0], [x[2]]) == 13x[1] + 2 + 2x[1]^2
end

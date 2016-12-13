@testset "Expectation" begin
    @polyvar x[1:3]
    p = x[3] - 2x[1]*x[2]^2 + 3x[3]*x[1] - 5x[1]^3
    v = [1,2,3]
    m = ζ(v, monomials(p), x)
    @test expectation(m, p) == p(v, x) == expectation(p, m)
    @test_throws ErrorException expectation(m, x[1] * x[2] * x[3])
    @test expectation(m, 0.5 * x[1] * x[2]^2) == 2.0
end

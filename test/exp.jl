@testset "Expectation" begin
    @polyvar x[1:3]
    p = x[3] - 2x[1]*x[2]^2 + 3x[3]*x[1] - 5x[1]^3
    v = (1,2,3)
    m = ζ(monomials(p), x=>v)
    @test expectation(m, p) == p(x=>v) == expectation(p, m)
    @test_throws ErrorException dot(x[1] * x[2] * x[3], m)
    @test dot(0.5 * x[1] * x[2]^2, m) == 2.0
    @test dot(m, x[1] * x[3]) == 3
end

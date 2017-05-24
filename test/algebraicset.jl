@testset "Algebraic set" begin
    @polyvar x y
    @test isa(FullSpace(), FullSpace)
    V = AlgebraicSet()
    addequality!(V, x*y - 1)
    @test_throws ArgumentError addinequality!(V, x*y)
    S = BasicSemialgebraicSet()
    addequality!(S, x)
    addinequality!(S, x^2*y)
    addequality!(S, 2*x^2*y)
    @test typeof(Int32(2)*x^2*y) == Term{true, Int32}
    addequality!(S, Int32(2)*x^2*y)
    addinequality!(S, 1.0x^2*y)
    addequality!(S, (6//3)*x^2*y+y)
    addinequality!(S, 1.5x+y)
    S2 = S ∩ V
    S3 = V ∩ S
    @test S2.p == S3.p == S.p
    @test S2.V.p == S3.V.p
    T = BasicSemialgebraicSet()
    addequality!(T, x*y^2 + 1)
    addinequality!(T, 1 - (x^2 + y^2))
    V2 = T.V ∩ V
    @test V2.p == [T.V.p; V.p]
    S4 = S ∩ T
    @test S4.p == [S.p; T.p]
    @test S4.V.p == [S.V.p; T.V.p]
end

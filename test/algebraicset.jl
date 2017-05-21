@testset "Algebraic set" begin
    @polyvar x y
    @test isa(FullSpace(), FullSpace)
    V = AlgebraicSet()
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
end

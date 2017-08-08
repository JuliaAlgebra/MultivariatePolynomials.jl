@testset "Algebraic set" begin
    Mod.@polyvar x y
    @test isa(FullSpace(), FullSpace)
    V = @set x * y == 1
    @test V isa AlgebraicSet{Int}
    @test_throws ArgumentError addinequality!(V, x*y)
    S = @set x == 0 && x^2*y >= 1
    addequality!(S, x^2 - y)
    addinequality!(S, x + y - 1)
    @test S isa BasicSemialgebraicSet{Int}
    @test Int32(2)*x^2*y isa AbstractTerm{Int32}
    S = (@set Int32(2)*x^2*y == 0 && 1.0x^2*y >= 0 && (6//3)*x^2*y == -y && 1.5x+y >= 0)
    S2 = S ∩ V
    S3 = V ∩ S
    @test S2.p == S3.p == S.p
    @test S2.V.p == S3.V.p
    T = (@set x*y^2 == -1 && x^2 + y^2 <= 1)
    V2 = @set T.V && V && x + y == 2.0
    @test V2 isa AlgebraicSet
    @test V2.p == [x + y - 2.0; T.V.p; V.p]
    S4 = @set S && T
    @test S4.p == [S.p; T.p]
    @test S4.V.p == [S.V.p; T.V.p]
end

@testset "NCMonomial" begin
    @ncpolyvar x
    @test_throws ArgumentError NCMonomial([x], [1,0])
    X = NCMonomial()
    @test nvars(X) == 0
    @test isempty(X.vars)
    @test isempty(X.z)
    X = NCMonomial(x)
    @test nvars(X) == 1
    @test X.vars[1] == x
    @test X.z[1] == 1
end
@testset "NCMonomialVector" begin
    @ncpolyvar x y
    @test_throws ArgumentError NCMonomialVector([x], [[1], [1,0]])
    X = NCMonomialVector()
    @test isempty(X.vars)
    @test isempty(X.Z)
    X = NCMonomialVector([x, y], [[1, 0], [0, 1]])
    @test X[1] == x
    @test X[2] == y
end

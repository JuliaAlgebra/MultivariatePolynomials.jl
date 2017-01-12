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
    X = NCMonomialVector([x,y], 2)
    @test X.vars == [x, y, x]
    @test X.Z == [[2, 0, 0], [1, 1, 0], [0, 2, 0], [0, 1, 1]]
    X0 = [x^3, x^2*y, x*y^2, x*y*x, y^3, y^2*x, y*x^2, y*x*y]
    X1 = monomials([x,y], 3)
    X2 = NCMonomialVector(X0)
    Z = [[3,0,0,0],[2,1,0,0],[1,2,0,0],[1,1,1,0],[0,3,0,0],[0,2,1,0],[0,1,2,0],[0,1,1,1]]
    @test typeof(X1) == Vector{NCMonomial}
    @test typeof(X2) == NCMonomialVector
    @test length(X1) == length(X2) == length(Z)
    @test issorted(X1, rev=true)
    @test issorted(X2, rev=true)
    @test issorted(Z, rev=true)
    @test X1 == X0
    for i in 1:length(Z)
        @test X1[i].vars == X2[i].vars == [x, y, x, y]
        @test X1[i].z == X2[i].z == Z[i]
    end
end

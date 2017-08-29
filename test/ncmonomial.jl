@testset "Non-commutative Monomial" begin
    Mod.@ncpolyvar x
    @test_throws ArgumentError Monomial{false}([x], [1,0])
    X = constantmonomial(x)
    @test nvariables(X) == 0
    @test isempty(X.vars)
    @test isempty(X.z)
    X = monomial(x)
    @test nvariables(X) == 1
    @test variables(X)[1] == x
    @test collect(exponents(X)) == [1]
    X = x^0
    Y = X * x^2
    @test nvariables(Y) == 1
    @test variables(Y)[1] == x
    @test collect(exponents(Y)) == [2]
    Y = x^2 * X
    @test nvariables(Y) == 1
    @test variables(Y)[1] == x
    @test collect(expoents(Y)) == [2]
end
@testset "Non-commutative MonomialVector" begin
    @ncpolyvar x y
    @test_throws ArgumentError MonomialVector{false}([x], [[1], [1,0]])
    X = MonomialVector{false}()
    @test isempty(X.vars)
    @test isempty(X.Z)
    X = MonomialVector{false}([x, y], [[1, 0], [0, 1]])
    @test X[1] == x
    @test X[2] == y
    X = MonomialVector([x, y], 2)
    @test X.vars == [x, y, x]
    @test X.Z == [[2, 0, 0], [1, 1, 0], [0, 2, 0], [0, 1, 1]]
    X0 = [x^3, x^2*y, x*y^2, x*y*x, y^3, y^2*x, y*x^2, y*x*y]
    X1 = monomials([x, y], 3)
    X2 = MonomialVector(X0)
    Z = [[3,0,0,0],[2,1,0,0],[1,2,0,0],[1,1,1,0],[0,3,0,0],[0,2,1,0],[0,1,2,0],[0,1,1,1]]
    @test typeof(X1) == Vector{Monomial{false}}
    @test typeof(X2) == MonomialVector{false}
    @test length(X1) == length(X2) == length(Z)
    @test issorted(X1, rev=true)
    @test issorted(X2, rev=true)
    @test issorted(Z, rev=true, lt=DynamicPolynomials.grlex)
    @test X1 == X0
    for i in 1:length(Z)
        @test X1[i].vars == X2[i].vars == [x, y, x, y]
        @test X1[i].z == X2[i].z == Z[i]
    end
end
@testset "PolyVar * Monomial" begin
    @ncpolyvar x y z
    m = y * Monomial([y, z, x, z], [0, 0, 2, 1])
    @test variables(m) == [y, z, x, z]
    @test m.z == [1, 0, 2, 1]
    m = x * Monomial([z, y, y, z], [0, 0, 2, 1])
    @test variables(m) == [z, x, y, y, z]
    @test m.z == [0, 1, 0, 2, 1]
    m = x * Monomial([y, z, y, z], [0, 0, 2, 1])
    @test variables(m) == [y, z, x, y, z]
    @test m.z == [0, 0, 1, 2, 1]
end
@testset "Monomial * PolyVar" begin
    @ncpolyvar x y z
    m = Monomial([x, z, x, y], [2, 1, 0, 0]) * y
    @test variables(m) == [x, z, x, y]
    @test m.z == [2, 1, 0, 1]
    m = Monomial([x, y, y, x], [2, 1, 0, 0]) * z
    @test variables(m) == [x, y, y, z, x]
    @test m.z == [2, 1, 0, 1, 0]
    m = Monomial([x, y, x, y], [2, 1, 0, 0]) * z
    @test variables(m) == [x, y, z, x, y]
    @test m.z == [2, 1, 1, 0, 0]
end

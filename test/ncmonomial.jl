@testset "Non-commutative Monomial" begin
    Mod.@ncpolyvar x
    X = constantmonomial(typeof(x))
    @test nvariables(X) == 0
    @test isempty(X.vars)
    @test isempty(X.z)
    X = constantmonomial(x)
    @test nvariables(X) == 1
    @test variables(X)[1] == x
    @test collect(exponents(X)) == [0]
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
    @test collect(exponents(Y)) == [2]
end
@testset "Non-commutative MonomialVector" begin
    Mod.@ncpolyvar x y
    X = emptymonovec(typeof(x))
    @test iszero(nvariables(X))
    @test isempty(variables(X))
    X = emptymonovec(x)
    @test nvariables(X) == 1
    @test variables(X)[1] == x
    X = monovec([y, x])
    @test X[1] == x
    @test X[2] == y
    X = monomials([x, y], 2)
    @test variables(X)[1] == x
    @test variables(X)[2] == y
    @test variables(X)[3] == x
    @test collect(exponents(X[1])) == [2, 0, 0]
    @test collect(exponents(X[2])) == [1, 1, 0]
    @test collect(exponents(X[3])) == [0, 2, 0]
    @test collect(exponents(X[4])) == [0, 1, 1]
    X0 = [x^3, x^2*y, x*y^2, x*y*x, y^3, y^2*x, y*x^2, y*x*y]
    X1 = monomials([x, y], 3)
    X2 = monovec(X0)
    Z = [[3,0,0,0],[2,1,0,0],[1,2,0,0],[1,1,1,0],[0,3,0,0],[0,2,1,0],[0,1,2,0],[0,1,1,1]]
    @test X1 isa AbstractVector{<:AbstractMonomial}
    @test X2 isa AbstractVector{<:AbstractMonomial}
    @test length(X1) == length(X2) == length(Z)
    @test issorted(X1, rev=true)
    @test issorted(X2, rev=true)
    @test X1 == X0
    for i in 1:length(Z)
        @test variables(X1[i])[1] == variables(X2[i])[1] == x
        @test variables(X1[i])[2] == variables(X2[i])[2] == y
        @test variables(X1[i])[3] == variables(X2[i])[3] == x
        @test variables(X1[i])[4] == variables(X2[i])[4] == y
        @test collect(exponents(X1[i])) == collect(exponents(X2[i])) == Z[i]
    end
end
#@testset "PolyVar * Monomial" begin
#    Mod.@ncpolyvar x y z
#    m = y * Monomial([y, z, x, z], [0, 0, 2, 1])
#    @test variables(m) == [y, z, x, z]
#    @test m.z == [1, 0, 2, 1]
#    m = x * Monomial([z, y, y, z], [0, 0, 2, 1])
#    @test variables(m) == [z, x, y, y, z]
#    @test m.z == [0, 1, 0, 2, 1]
#    m = x * Monomial([y, z, y, z], [0, 0, 2, 1])
#    @test variables(m) == [y, z, x, y, z]
#    @test m.z == [0, 0, 1, 2, 1]
#end
#@testset "Monomial * PolyVar" begin
#    @ncpolyvar x y z
#    m = Monomial([x, z, x, y], [2, 1, 0, 0]) * y
#    @test variables(m) == [x, z, x, y]
#    @test m.z == [2, 1, 0, 1]
#    m = Monomial([x, y, y, x], [2, 1, 0, 0]) * z
#    @test variables(m) == [x, y, y, z, x]
#    @test m.z == [2, 1, 0, 1, 0]
#    m = Monomial([x, y, x, y], [2, 1, 0, 0]) * z
#    @test variables(m) == [x, y, z, x, y]
#    @test m.z == [2, 1, 1, 0, 0]
#end

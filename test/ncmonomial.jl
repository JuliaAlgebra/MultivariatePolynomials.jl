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
    # `variables` may return `[x, y, x, y]` if the polynomial has, e.g., the monomial `x * y * x * y`.
    X = monomials([y, x, y, x, y], 0:2)
    Y = monomials([x, y], 2)
    for M in [X, Y]
        @test variables(X)[1] == x
        @test variables(X)[2] == y
        @test variables(X)[3] == x
        @test collect(exponents(X[1])) == [2, 0, 0]
        @test collect(exponents(X[2])) == [1, 1, 0]
        @test collect(exponents(X[3])) == [0, 2, 0]
        @test collect(exponents(X[4])) == [0, 1, 1]
    end
    @test collect(exponents(X[5])) == [1, 0, 0]
    @test collect(exponents(X[6])) == [0, 1, 0]
    @test collect(exponents(X[7])) == [0, 0, 0]
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

@testset "Non-commutative Monomial" begin
    Mod.@ncpolyvar x
    X_ = constantmonomial(typeof(x))
    @test nvariables(X_) == 0
    @test isempty(X_.vars)
    @test isempty(X_.z)
    X0 = constantmonomial(x)
    for Y in [x^0, X0]
        @test nvariables(Y) == 1
        @test variables(Y)[1] == x
        @test collect(exponents(Y)) == [0]
    end
    for X in [X_, X0]
        for Y in [monomial(x), x^1, x * X, X * x]
            @test nvariables(Y) == 1
            @test variables(Y)[1] == x
            @test collect(exponents(Y)) == [1]
        end
    end
    for X in [X_, X0]
        for Y in [x^2, X * x^2, x^2 * X]
            @test nvariables(Y) == 1
            @test variables(Y)[1] == x
            @test collect(exponents(Y)) == [2]
        end
    end
    @testset "Issue #148" begin
        Mod.@ncpolyvar x y z
        m = x * y^4 * x^2
        @test 3 == @inferred degree(m, x)
        @test 4 == @inferred degree(m, y)
        @test 0 == @inferred degree(m, z)
    end
    @testset "Issue #71 of DynamicPolynomials" begin
        Mod.@ncpolyvar x y
        @test x^0 * y == y * x^0
    end
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
    @test X[1] == y
    @test X[2] == x
    # `variables` may return `[x, y, x, y]` if the polynomial has, e.g., the monomial `x * y * x * y`.
    X = monomials([y, x, y, x, y], 0:2)
    Y = monomials([x, y], 2)
    for M in [X, Y]
        @test variables(X)[1] == x
        @test variables(X)[2] == y
        @test variables(X)[3] == x
        @test collect(exponents(X[7])) == [2, 0, 0]
        @test collect(exponents(X[6])) == [1, 1, 0]
        @test collect(exponents(X[5])) == [0, 2, 0]
        @test collect(exponents(X[4])) == [0, 1, 1]
    end
    @test collect(exponents(X[3])) == [1, 0, 0]
    @test collect(exponents(X[2])) == [0, 1, 0]
    @test collect(exponents(X[1])) == [0, 0, 0]
    X0 = [x^3, x^2*y, x*y^2, x*y*x, y^3, y^2*x, y*x^2, y*x*y]
    X1 = monomials([x, y], 3)
    X2 = monovec(X0)
    Z = reverse([[3,0,0,0],[2,1,0,0],[1,2,0,0],[1,1,1,0],[0,3,0,0],[0,2,1,0],[0,1,2,0],[0,1,1,1]])
    @test X1 isa AbstractVector{<:AbstractMonomial}
    @test X2 isa AbstractVector{<:AbstractMonomial}
    @test length(X1) == length(X2) == length(Z)
    @test issorted(X1)
    @test issorted(X2)
    @test X1 == X0
    for i in eachindex(Z)
        @test variables(X1[i])[1] == variables(X2[i])[1] == x
        @test variables(X1[i])[2] == variables(X2[i])[2] == y
        @test variables(X1[i])[3] == variables(X2[i])[3] == x
        @test variables(X1[i])[4] == variables(X2[i])[4] == y
        @test collect(exponents(X1[i])) == collect(exponents(X2[i])) == Z[i]
    end
end

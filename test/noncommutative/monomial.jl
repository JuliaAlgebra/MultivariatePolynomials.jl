if !@isdefined(Mod)
    import DynamicPolynomials as Mod
    include("../utils.jl")
end

@testset "Non-commutative Monomial" begin
    Mod.@ncpolyvar x
    X_ = constant_monomial(typeof(x))
    @test nvariables(X_) == 0
    @test isempty(variables(X_))
    @test isempty(exponents(X_))
    X0 = constant_monomial(x)
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
    v = Mod.@ncpolyvar x y
    @test_throws ArgumentError monomials(v, -1)
    @test_throws ArgumentError monomials(v, -1:1)
    @test_throws ArgumentError monomials(v, [1, -1])
    X = empty_monomial_vector(typeof(x))
    @test iszero(nvariables(X))
    @test isempty(variables(X))
    X = empty_monomial_vector(x)
    # An ordinary empty vector has no variables stored alongside its entries.
    @test iszero(nvariables(X))
    @test isempty(variables(X))
    X = monomial_vector([y, x])
    @test X[1] == y
    @test X[2] == x
    # `variables` may return `[x, y, x, y]` if the polynomial has, e.g., the monomial `x * y * x * y`.
    X = monomials([y, x, y, x, y], 0:2)
    Y = monomials([x, y], 2)
    for M in [X, Y]
        @test variables(M) == [x, y, x]
        @test collect(exponents(M[end])) == [2, 0, 0]
        @test collect(exponents(M[end-1])) == [1, 1, 0]
        @test collect(exponents(M[end-2])) == [0, 2, 0]
        @test collect(exponents(M[end-3])) == [0, 1, 1]
    end
    @test collect(exponents(X[3])) == [1, 0, 0]
    @test collect(exponents(X[2])) == [0, 1, 0]
    @test collect(exponents(X[1])) == [0, 0, 0]
    X0 = [x^3, x^2 * y, x * y^2, x * y * x, y^3, y^2 * x, y * x^2, y * x * y]
    original = copy(X0)
    X1 = monomials([x, y], 3)
    X2 = monomial_vector(X0)
    Z = reverse([
        [3, 0, 0, 0],
        [2, 1, 0, 0],
        [1, 2, 0, 0],
        [1, 1, 1, 0],
        [0, 3, 0, 0],
        [0, 2, 1, 0],
        [0, 1, 2, 0],
        [0, 1, 1, 1],
    ])
    @test X1 isa AbstractVector{<:AbstractMonomial}
    @test X2 isa AbstractVector{<:AbstractMonomial}
    @test length(X1) == length(X2) == length(Z)
    @test issorted(X1)
    @test issorted(X2)
    @test X0 == original
    @test X1 == X2 == reverse(X0)
    vars = [x, y, x, y]
    @test variables(X1) == variables(X2) == vars
    for i in eachindex(Z)
        @test variables(X1[i]) == vars
        # Sorting keeps each input monomial's variable sequence. Align it
        # explicitly to compare with the padded enumeration representation.
        @test collect(exponents(X1[i])) == exponents(X2[i], vars) == Z[i]
    end
end
@testset "is_commutative" begin
    function _is_comm_test(p, expected)
        @test (@inferred MP.is_commutative(p)) === expected
        @test (@inferred MP.is_commutative(typeof(p))) === expected
    end
    Mod.@ncpolyvar x y[1:2]
    for p in (x, y, x^2, 2x^2, sum(y))
        _is_comm_test(p, false)
    end
    Mod.@polyvar x y[1:2]
    for p in (x, y, x^2, 2x^2, sum(y))
        _is_comm_test(p, true)
    end
end

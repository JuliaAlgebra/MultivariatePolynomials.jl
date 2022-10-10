@testset "Monomial Vector" begin
    Mod.@polyvar x y
    @test x > y
    @test x^2 > y^2
    X = [y^2, x*y, x^2]
    @test isempty(@inferred monomials((x, y), 1:0))
    for (i, m) in enumerate(monomials((x, y), 2))
        @test m == X[i]
    end
    for (i, m) in enumerate(monomials((x, y), 2:2))
        @test m == X[i]
    end
    X = [y^2, x^2]
    for (i, m) in enumerate(monomials((x, y), 2, m -> m != x*y))
        @test m == X[i]
    end
    @test (@inferred monovectype([1, x])) <: AbstractArray{<:AbstractMonomial}
    @test (@inferred monovectype([x])) <: AbstractArray{<:AbstractMonomial}
    @test (@inferred monovec([1, x])) isa monovectype([1, x])
    @test (@inferred monovec([x])) isa monovectype([x])
    @test (@inferred monovec([1, 2], [1, x]))[2] isa AbstractArray{<:AbstractMonomial}
    @test (@inferred monovec([1], [x]))[2] isa AbstractArray{<:AbstractMonomial}
    @test length(monovec([y, x])) == 2
    X = monovec([x, 1, x*y])
    @test X == collect(X)
    @test nvariables(X) == 2
    @test variables(X)[1] == x
    @test variables(X)[2] == y
    @test X[2:3][1] == x
    @test X[2:3][2] == x * y
    @test monovec(X[[3, 2]])[1] == x
    @test monovec(X[[3, 2]])[2] == x * y
    # Documentation examples
    @test monovec([x*y, x, x*y, x^2*y, x]) == [x, x*y, x^2*y]
    @test monovectype([x*y, x, 1, x^2*y, x]) <: AbstractVector{typeof(x*y)}
    @test monovectype([x*y, x, x*y, x^2*y, x]) <: AbstractVector
    σ, smv = sortmonovec([x*y, x, x*y, x^2*y, x])
    @test smv == [x, x*y, x^2*y]
    @test σ[3] == 4
    @test σ[2] in (1, 3)
    @test σ[1] in (2, 5)
    @test mergemonovec([[x*y, x, x*y], [x^2*y, x]]) == [x, x*y, x^2*y]
    @test_throws ArgumentError monovec([1, 2], [x^2])
    σ, X = sortmonovec((x, y))
    @test σ == [2, 1]
    @test X == [y, x]
    @test monomialtype([x, y]) <: AbstractMonomial
    @test monomialtype([x^2, 1]) <: AbstractMonomial
    @test monomialtype([x*y, x+y]) <: AbstractMonomial

    @test monovec([x, x^2]) != monovec([x*y, x^2*y])
    @test monomials(x, 1:3) == monovec([x^3, x, x^2])
    @test monomials((x, y), 2) != monomials((x, y), 1)

    # See https://github.com/JuliaAlgebra/DynamicPolynomials.jl/issues/111
    @testset "Indexing with vector of boolean" begin
        vars = Mod.@polyvar x y
        X = monomials(vars, 2)
        @test X[[true, false, true]] == monovec([x^2, y^2])
        X = monomials(vars, 0:1)
        @test filter(mono -> degree(mono) == 1, X) == monovec([x, y])
        @test filter(mono -> degree(mono) == 0, X) == monovec([x^0])
    end
end

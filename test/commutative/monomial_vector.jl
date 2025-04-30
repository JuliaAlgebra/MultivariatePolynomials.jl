@testset "Monomial Vector" begin
    Mod.@polyvar x y
    @test x > y
    @test x^2 > y^2
    X = [y^2, x * y, x^2]
    @test isempty(@inferred monomials((x, y), 1:0))
    monoss = [monomials((x, y), 2), monomials((x, y), 2:2)]
    # TypedPolynomials is allowed to error on `monomials((y, x), 2)`
    # because it would make `monomials` type unstable to sort the tuple of variables
    # of different types
    if typeof(x) === typeof(y)
        push!(monoss, monomials((y, x), 2))
        push!(monoss, monomials((y, x), 2:2))
    end
    for monos in monoss
        for (i, m) in enumerate(monos)
            @test m == X[i]
        end
    end
    X = [y^2, x^2]
    for (i, m) in enumerate(monomials((x, y), 2, m -> m != x * y))
        @test m == X[i]
    end
    @test (@inferred monomial_vector_type([1, x])) <:
          AbstractArray{<:AbstractMonomial}
    @test (@inferred monomial_vector_type([x])) <:
          AbstractArray{<:AbstractMonomial}
    @test (@inferred monomial_vector([1, x])) isa monomial_vector_type([1, x])
    @test (@inferred monomial_vector([x])) isa monomial_vector_type([x])
    @test (@inferred monomial_vector([1, 2], [1, x]))[2] isa
          AbstractArray{<:AbstractMonomial}
    @test (@inferred monomial_vector([1], [x]))[2] isa
          AbstractArray{<:AbstractMonomial}
    @test length(monomial_vector([y, x])) == 2
    X = monomial_vector([x, 1, x * y])
    @test X == collect(X)
    @test nvariables(X) == 2
    @test variables(X)[1] == x
    @test variables(X)[2] == y
    @test X[2:3][1] == x
    @test X[2:3][2] == x * y
    @test monomial_vector(X[[3, 2]])[1] == x
    @test monomial_vector(X[[3, 2]])[2] == x * y
    # Documentation examples
    @test monomial_vector([x * y, x, x * y, x^2 * y, x]) == [x, x * y, x^2 * y]
    @test monomial_vector_type([x * y, x, 1, x^2 * y, x]) <:
          AbstractVector{typeof(x * y)}
    @test monomial_vector_type([x * y, x, x * y, x^2 * y, x]) <: AbstractVector
    σ, smv = sort_monomial_vector([x * y, x, x * y, x^2 * y, x])
    @test smv == [x, x * y, x^2 * y]
    @test σ[3] == 4
    @test σ[2] in (1, 3)
    @test σ[1] in (2, 5)
    @test merge_monomial_vectors([[x * y, x, x * y], [x^2 * y, x]]) ==
          [x, x * y, x^2 * y]
    @test_throws ArgumentError monomial_vector([1, 2], [x^2])
    σ, X = sort_monomial_vector((x, y))
    @test σ == [2, 1]
    @test X == [y, x]
    @test monomial_type([x, y]) <: AbstractMonomial
    @test monomial_type([x^2, 1]) <: AbstractMonomial
    @test monomial_type([x * y, x + y]) <: AbstractMonomial

    @test monomial_vector([x, x^2]) != monomial_vector([x * y, x^2 * y])
    @test monomials(x, 1:3) == monomial_vector([x^3, x, x^2])
    @test monomials((x, y), 2) != monomials((x, y), 1)

    # See https://github.com/JuliaAlgebra/DynamicPolynomials.jl/issues/111
    @testset "Indexing with vector of boolean" begin
        vars = Mod.@polyvar x y
        X = monomials(vars, 2)
        @test X[[true, false, true]] == monomial_vector([x^2, y^2])
        X = monomials(vars, 0:1)
        @test filter(mono -> degree(mono) == 1, X) == monomial_vector([x, y])
        @test filter(mono -> degree(mono) == 0, X) == monomial_vector([x^0])
    end

    @testset "monomials" begin
        Mod.@polyvar v[1:3]
        @test monomials(v, 0:3) == [
            v[1]^0,
            v[3],
            v[2],
            v[1],
            v[3]^2,
            v[2] * v[3],
            v[2]^2,
            v[1] * v[3],
            v[1] * v[2],
            v[1]^2,
            v[3]^3,
            v[2] * v[3]^2,
            v[2]^2 * v[3],
            v[2]^3,
            v[1] * v[3]^2,
            v[1] * v[2] * v[3],
            v[1] * v[2]^2,
            v[1]^2 * v[3],
            v[1]^2 * v[2],
            v[1]^3,
        ]
    end

    @testset "Negative degree" begin
        Mod.@polyvar x[1:3]
        @test_throws ArgumentError monomials(x, -1)
        @test_throws ArgumentError monomials(x, -1:1)
        @test_throws ArgumentError monomials(x, [1, -1])
    end
end

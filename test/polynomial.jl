using Test
using LinearAlgebra
import MutableArithmetics
const MA = MutableArithmetics
using MultivariatePolynomials
const MP = MultivariatePolynomials
@testset "Polynomial" begin
    Mod.@polyvar x

    @test terms(polynomial([1, x^2, x, 2x^2])) == [1, x, 3x^2]
    @test terms(polynomial([x, 3x^4, 2], MP.UniqState())) == [2, x, 3x^4]
    @test terms(polynomial([-2, 2, x, x^2, -2x^2, x^2, x^3, 2x^3], MP.SortedState())) == [x, 3x^3]

    @test polynomial(1 + x) == 1 + x
    @test leadingterm(1 + x) == x
    @test leadingterm(x - x) == 0
    @test leadingmonomial(x - x) == 1
    @test leadingcoefficient(x - x) == 0
    @test one(1 + x) == one(1.0 + x) == 1
    @test zero(1 + x) == zero(1.0 + x) == 0
    @test 1 != 1 + x
    @test 2x == x + x^2 + x - x^2
    @test x + x^2 - x^2 - x != 2x
    @test x^2 + x != x^2 + x + 1
    #@inferred one(1 + x)
    @inferred zero(1 + x)
    #@inferred one(1.0 + x)
    @inferred zero(1.0 + x)

    @test polynomial(1:2, [x, x^2]) == x + 2x^2
    @test polynomial(1:2, monomials(x, 1:2)) == x + 2x^2

    @test terms(polynomial(1 + x + x^2 - x + x^2)) == [1, 2x^2]
    @test terms(CustomPoly(1 + x + x^2 - x + x^2)) == [1, 2x^2]

    @test (1.0 + x) * x == x^2 + x
    @test constantterm(1, x) * (1 - x) == 1 - x
    @test promote_type(typeof(1-x), typeof(x)) <: AbstractPolynomial{Int}
    @test x != 1 - x

    @test term(x + x^2 - x) isa AbstractTerm
    @test term(x + x^2 - x) == x^2
    @test term(x - x) isa AbstractTerm
    @test iszero(term(x - x))
    @test_throws InexactError term(x + x^2)

    Mod.@polyvar y

    p = 3x^2*y^4 + 2x
    @test +(p) === p
    @test *(p) === p
    @test terms(p)[1] == 2x
    @test terms(p)[end] == 3x^2*y^4
    typetests(p)
    typetests([p, x + y])
    @test (@inferred polynomial(p)) isa AbstractPolynomial{Int}
    @test (@inferred polynomial(p, Float64)) isa AbstractPolynomial{Float64}

    @test coefficient(2x + 4y^2 + 3, y^2) == 4
    @test coefficient(2x + 4y^2 + 3, x^2) == 0

    Mod.@polyvar a b
    @test coefficient((2a + b)x^2 + 2y^2 + 3, x^2, (x,y)) == 2a + b
    @test coefficient((2a + b)x^2 + 2y^2 + 3x^2*y, x^2, (x,)) == 2a + b + 3y

    @test (@inferred 2x^2*y + 0.0x*y) == 2x^2*y
    @test (@inferred 0.0x^2*y + 3x*y) == 3x*y

    @test iszero(((x + x) - 2x) * (x * (x ^ 2 + y ^ 2)))

    @test Tuple(variables([x + 1, y^2])) == (x, y)
    @test Tuple(variables([y^2, x + 1])) == (x, y)

    @test maxdegree(x^2 - x^2) == 0
    @test maxdegree(x^2 - x^2, x) == 0
    @test maxdegree(x^2 - x^2, y) == 0
    @test mindegree(x^2 - x^2) == 0
    @test mindegree(x^2 - x^2, x) == 0
    @test mindegree(x^2 - x^2, y) == 0
    @test extdegree(x^2 - x^2) == (0, 0)
    @test extdegree(x^2 - x^2, x) == (0, 0)
    @test extdegree(x^2 - x^2, y) == (0, 0)
    @test maxdegree(x*y + 2 + x^2*y + x + y) == 3
    @test maxdegree(x*y + 2 + x^2*y + x + y, x) == 2
    @test maxdegree(x*y + 2 + x^2*y + x + y, y) == 1
    @test mindegree(x*y + 2 + x^2*y + x + y) == 0
    @test mindegree(x*y + 2 + x^2*y + x + y, x) == 0
    @test mindegree(x*y + 2 + x^2*y + x + y, y) == 0
    @test extdegree(x*y + 2 + x^2*y + x + y) == (0, 3)
    @test extdegree(x*y + 2 + x^2*y + x + y, x) == (0, 2)
    @test extdegree(x*y + 2 + x^2*y + x + y, y) == (0, 1)
    @test extdegree(x*y + x^2*y, x) == (1, 2)
    @test extdegree(x*y + x^2*y, y) == (1, 1)
    @test leadingterm(x*y + 2 + x^2*y + x + y) == x^2*y
    @test nvariables(x + y - x) == 2
    @test nvariables(x + x^2) == 1

    @test collect(coefficients(x*y + 2 + 3x^2*y + 4x + 6y, [x, x*y^2, x*y, x^2*y, y, x^3])) == [4, 0, 1, 3, 6, 0]

    # Doc examples
    @test collect(coefficients(4x^2*y + x*y + 2x)) == [2, 1, 4]
    @test collect(coefficients(4x^2*y + x*y + 2x + 3, [x, 1, x*y, y])) == [2, 3, 1, 0]

    for p in [polynomial([4, 9], [x, x*x]), polynomial([9, 4], [x*x, x])]
        @test collect(coefficients(p)) == [4, 9]
        @test monomials(p)[1] == x
        @test monomials(p)[2] == x^2
        @test p == dot([4, 9], [x, x*x])
    end

    @inferred polynomial(i -> float(i), [x, x*x])
    @inferred polynomial(i -> 3 - float(i), monovec([x*x, x]))
    for p in [polynomial(i -> float(i), [x, x*x]),
              polynomial(i -> 1.0, [x*x, x, x*x]),
              polynomial(i -> float(i), monovec([x*x, x]))]
        @test collect(coefficients(p)) == [1.0, 2.0]
        @test collect(monomials(p)) == monovec([x, x^2])
    end

    @test (x + y)' == x + y
    @test transpose(x + y) == x + y
    @test transpose([1 2; 3 4] * x) == [1 3; 2 4] * x

    @test removemonomials(4x^2*y + x*y + 2x, [x*y]) == 4x^2*y + 2x

    @test_throws InexactError push!([1], x+1)

    @test polynomial([1 2; 3 4], [x^2, y]) == x^4 + 5x^2*y + 4y^2
    @test polynomial([1 2; 3 4], [x^2, y], Float64) isa AbstractPolynomial{Float64}
    @test polynomial([1 2; 3 4], [y, x^2]) == y^2 + 5x^2*y + 4x^4
    @test polynomial([1 2; 3 4], [y, x^2], Float64) isa AbstractPolynomial{Float64}
    @test polynomial([1 2; 3 4], monovec([y, x^2])) == 4x^4 + 5x^2*y + y^2
    @test polynomial([1 2; 3 4], monovec([y, x^2]), Float64) isa AbstractPolynomial{Float64}

    @test (@inferred round(2.6x + 1.001x^2)) == 3x + 1x^2
    @test (@inferred round(3.1x*y)) == 3x*y
    @test (@inferred round(2.613x + 1.1051x^2, digits=2)) ≈ 2.61x + 1.11x^2
    @test (@inferred round(3.145x*y, digits=1)) ≈ 3.1x*y

    @testset "Graded Lex Order" begin
        Mod.@polyvar x y z
        p = 3*y^2 + 2*y*x
        @test collect(coefficients(p)) == [3, 2]
        @test collect(monomials(p)) == monovec([x*y, y^2])
        # Examples from p. 59 of the 4th edition of "Ideals, Varieties, and Algorithms" of Cox, Little and O'Shea
        f = 4*x*y^2*z + 4*z^2 - 5*x^3 + 7*x^2*z^2
        @test collect(coefficients(f)) == [4, -5, 4, 7]
        @test collect(monomials(f)) == monovec([x^2*z^2, x*y^2*z, x^3, z^2])
        @test ordering(f) === GradedLex()
    end

    @testset "Convertion" begin
        Mod.@polyvar x y z
        p = 2.5x + 1 - 2.5x
        @test convert(Int, p) == 1
        @test convert(typeof(p), p) === p
        @test convert(Union{Nothing, typeof(p)}, p) === p
        a = 2y
        q = polynomial([a, z, -a], [x, 1, x])
        @test convert_to_constant(q) == z
        q = polynomial([a, z], [x, 1])
        @test_throws InexactError convert_to_constant(q)
        alloc_test(() -> convert(typeof(p), p), 0)
    end

    @testset "Vector" begin
        Mod.@polyvar x y
        v = [x - 1, y + 1]
        @test nvariables(v) == 2
        @test variables(v)[1] == x
        @test variables(v)[2] == y
    end

    @testset "Effective variables" begin
        Mod.@polyvar x y z
        T = variable_union_type(x)
        @test x isa T
        @test y isa T
        @test z isa T
        @test T[x] == @inferred effective_variables(x + y - y)
        @test T[x, y] == @inferred effective_variables(z + x + y - z)
        @test T[y, z] == @inferred effective_variables(z + 0 * x + y)
        @test T[z] == @inferred effective_variables(z + 0 * x + y^0)
    end

    @testset "Complex" begin
        # See https://github.com/JuliaAlgebra/MultivariatePolynomials.jl/issues/128
        Mod.@polyvar x
        p = im * x + 2im * x^2
        @test p * p == -x^2 - 4x^4 - 4x^3
    end

    @testset "$f is a mutable copy, see issue DynamicPolynomials#62" for f in [zero, one]
        p = 2x + 1
        q = f(p)
        q = MA.add!!(q, 2y)
        @test p == 2x + 1
    end

    @testset "mapcoefficients $nz" for nz in [false, true]
        p = 2x + 1
        @test mapcoefficients(x -> x / 2, p, nonzero = nz) == 1.0x + 0.5
        @test mapcoefficients(x -> x / 2, 3x, nonzero = nz) == 1.5x
        @test p === mapcoefficients!(x -> x + 1, p, nonzero = nz)
        @test p == 3x + 2
        q = zero(p)
        @test q === mapcoefficients_to!(q, x -> 2x, p, nonzero = nz)
        @test q == 6x + 4
        @test q === mapcoefficients_to!(q, x -> 2x, 3x, nonzero = nz)
        @test q == 6x
        @test q === mapcoefficients_to!(q, x -> 2x, x, nonzero = nz)
        @test q == 2x
        @test q === mapcoefficients_to!(q, x -> 2x, x^2, nonzero = nz)
        @test q == 2x^2
    end
end

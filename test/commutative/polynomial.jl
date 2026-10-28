if !@isdefined(Mod)
    import DynamicPolynomials as Mod
    include("../utils.jl")
end

using Test
using LinearAlgebra
import MutableArithmetics as MA
import StarAlgebras as SA
using MultivariatePolynomials
const MP = MultivariatePolynomials

@testset "Polynomial coefficient and monomial lengths" begin
    Mod.@polyvar x
    @test_throws ArgumentError polynomial([1, 2], [x])
    @test_throws ArgumentError polynomial([1], [x, x^2])
    @test_throws ArgumentError polynomial(Int[], [x])
end

@testset "Accumulate monomials in a prepared polynomial" begin
    Mod.@polyvar x y
    p = x + y + 1
    for op in (+, -), m in (x, y^2)
        output = copy(p)
        expected = op(p, m)
        @test (@inferred MA.operate!!(op, output, m)) === output
        @test output == expected
        @test p == x + y + 1
    end
end

@testset "Mutable operations promote insufficient destination bases" begin
    Mod.@polyvar x y
    for T in (Int, BigInt)
        original = T(2) * x + T(1)
        for (op, args) in (
            (+, (y,)),
            (-, (T(2) * y,)),
            (*, (T(2) * y + T(1),)),
            (MA.add_mul, (T(2) * y, y + T(1))),
            (MA.sub_mul, (y + T(1), T(2) * y)),
            (MA.add_mul, (transpose([T(2) * y]), [y + T(1)])),
        )
            expected = op(original, args...)
            saved_args = MA.copy_if_mutable.(args)
            mutable =
                MA.mutability(
                    typeof(original),
                    op,
                    typeof(original),
                    typeof.(args)...,
                ) isa MA.IsMutable
            output = MA.mutable_copy(original)
            if mutable
                @test (@inferred MA.operate!(op, output, args...)) === output
                @test output == expected
            else
                @test_throws ArgumentError MA.operate!(op, output, args...)
                @test output == original
            end
            buffer = MA.buffer_for(op, typeof(original), typeof.(args)...)
            for call in (
                output -> MA.operate!!(op, output, args...),
                output -> MA.operate_to!!(output, op, original, args...),
                output -> MA.buffered_operate!!(buffer, op, output, args...),
                output -> MA.buffered_operate_to!!(
                    buffer,
                    output,
                    op,
                    original,
                    args...,
                ),
            )
                output = MA.mutable_copy(original)
                alg = parent(output)
                @test MA.mutability(output, op, output, args...) ===
                      MA.mutability(
                    typeof(output),
                    op,
                    typeof(output),
                    typeof.(args)...,
                )
                result = @inferred call(output)
                @test result == expected
                @test (result === output) == mutable
                @test output == (mutable ? expected : original)
                @test parent(original) === alg
                @test parent(result) == parent(expected)
                if !mutable
                    @test parent(output) === alg
                end
                @test args == saved_args
            end

            output = original + zero(T(1) * x + T(1) * y)
            @test MA.mutability(output, op, output, args...) isa MA.IsMutable
            @test (@inferred MA.operate!!(op, output, args...)) === output
            @test output == expected
            @test args == saved_args
        end
        @test MA.@rewrite((original + 1) + y) == original + 1 + y
        @test MA.@rewrite(original + 2y * (y + 1)) == original + 2y * (y + 1)
        @test original == T(2) * x + T(1)
    end
end

@testset "Polynomial" begin
    Mod.@polyvar x

    @test terms(polynomial([1, x^2, x, 2x^2])) == [1, x, 3x^2]
    @test terms(polynomial([x, 3x^4, 2], MP.UniqState())) == [2, x, 3x^4]
    @test terms(
        polynomial([-2, 2, x, x^2, -2x^2, x^2, x^3, 2x^3], MP.SortedState()),
    ) == [x, 3x^3]

    @test polynomial(1 + x) == 1 + x
    @test leading_term(1 + x) == x
    @test leading_term(x - x) == 0
    @test leading_monomial(x - x) == 1
    @test leading_coefficient(x - x) == 0
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

    @testset "Coefficient vector types" begin
        for T in (Int, BigInt, Real, Number),
            (c, m, expected) in (
                (T[1, 2], [x, x^2], x + 2x^2),
                (T[1, -1], [x, x], 0),
                (T[], typeof(x)[], 0),
            )

            p = @inferred polynomial(c, m)
            @test p == expected
            @test coefficient_type(p) === T
        end
        @test coefficient_type(polynomial(Number[1, 0.5], [x, x^2])) === Number
    end

    @test terms(polynomial(1 + x + x^2 - x + x^2)) == [1, 2x^2]

    @testset "Scalar coefficient promotion" begin
        for (c, a) in ((2, 0.5), (2.0, 5.0im), ([1 2; 3 4], [2 3; 4 5])),
            (p, k) in
            ((term(c, x), zero(c)), (polynomial([c, c], [one(x), x]), c))

            original = copy(p)
            for (q, cx, c0) in (
                ((@inferred p + a), c, k + a),
                ((@inferred a + p), c, a + k),
                ((@inferred p - a), c, k - a),
                ((@inferred a - p), -c, a - k),
            )
                @test coefficient_type(q) === typeof(c + a)
                @test q == polynomial([c0, cx], [one(x), x])
            end
            @test p == original
        end
    end

    @test (1.0 + x) * x == x^2 + x
    @test constant_term(1, x) * (1 - x) == 1 - x
    @test promote_type(typeof(1 - x), typeof(x)) <: AbstractPolynomial{Int}
    @test x != 1 - x

    @test term(x + x^2 - x) isa AbstractTerm
    @test term(x + x^2 - x) == x^2
    @test term(x - x) isa AbstractTerm
    @test iszero(term(x - x))
    @test_throws InexactError term(x + x^2)

    Mod.@polyvar y

    p = 3x^2 * y^4 + 2x
    @test +(p) === p
    @test *(p) === p
    @test terms(p)[1] == 2x
    @test terms(p)[end] == 3x^2 * y^4
    typetests(p, Int)
    typetests([p, x + y], Int)
    @test (@inferred polynomial(p)) isa AbstractPolynomial{Int}
    @test (@inferred polynomial(p, Float64)) isa AbstractPolynomial{Float64}

    @test coefficient(2x + 4y^2 + 3, y^2) == 4
    @test coefficient(2x + 4y^2 + 3, x^2) == 0

    Mod.@polyvar a b
    @test coefficient((2a + b)x^2 + 2y^2 + 3, x^2, (x, y)) == 2a + b
    @test coefficient((2a + b)x^2 + 2y^2 + 3x^2 * y, x^2, (x,)) == 2a + b + 3y

    @test (@inferred 2x^2 * y + 0.0x * y) == 2x^2 * y
    @test (@inferred 0.0x^2 * y + 3x * y) == 3x * y

    @test iszero(((x + x) - 2x) * (x * (x^2 + y^2)))

    @test Tuple(variables([x + 1, y^2])) == (x, y)
    @test Tuple(variables([y^2, x + 1])) == (x, y)

    p = x^2 - x^2
    @test maxdegree(p) == 0
    @test maxdegree(p, x) == 0
    @test maxdegree(p, y) == 0
    @test mindegree(p) == 0
    @test mindegree(p, x) == 0
    @test mindegree(p, y) == 0
    @test extdegree(p) == (0, 0)
    @test extdegree(p, x) == (0, 0)
    @test extdegree(p, y) == (0, 0)
    q = x * y + 2 + x^2 * y + x + y
    @test maxdegree(q) == 3
    @test maxdegree(q, x) == 2
    @test maxdegree(q, y) == 1
    @test mindegree(q) == 0
    @test mindegree(q, x) == 0
    @test mindegree(q, y) == 0
    @test extdegree(q) == (0, 3)
    @test extdegree(q, x) == (0, 2)
    @test extdegree(q, y) == (0, 1)
    r = x * y + x^2 * y
    @test extdegree(r, x) == (1, 2)
    @test extdegree(r, y) == (1, 1)
    @test leading_term(q) == x^2 * y
    @test nvariables(x + y - x) == 2
    @test nvariables(x + x^2) == 1
    @test mindegree([p, q, r]) == 0
    @test mindegree([p, q, r], x) == 0
    @test mindegree([p, q, r], y) == 0
    @test mindegree([p, q], y) == 0
    @test mindegree([q], y) == 0
    @test mindegree([r], y) == 1
    @test maxdegree([p, q, r]) == 3
    @test maxdegree([p, q, r], x) == 2
    @test maxdegree([p, q, r], y) == 1
    @test extdegree([p, q, r]) == (0, 3)
    @test extdegree([p, q, r], x) == (0, 2)
    @test extdegree([p, q, r], y) == (0, 1)

    @test collect(
        coefficients(
            x * y + 2 + 3x^2 * y + 4x + 6y,
            [x, x * y^2, x * y, x^2 * y, y, x^3],
        ),
    ) == [4, 0, 1, 3, 6, 0]

    # Doc examples
    @test collect(coefficients(4x^2 * y + x * y + 2x)) == [2, 1, 4]
    @test collect(coefficients(4x^2 * y + x * y + 2x + 3, [x, 1, x * y, y])) ==
          [2, 3, 1, 0]
    @test monomials(4x^2 * y + x * y + 2x + 3)[1:1] ==
          [constant_monomial(x * y)]

    for p in [polynomial([4, 9], [x, x * x]), polynomial([9, 4], [x * x, x])]
        @test collect(coefficients(p)) == [4, 9]
        @test monomials(p)[1] == x
        @test monomials(p)[2] == x^2
        @test monomials(p)[1:2][1] == x
        @test monomials(p)[1:2][2] == x^2
        @test p == dot([4, 9], [x, x * x])
    end

    @inferred polynomial(i -> float(i), [x, x * x])
    @inferred polynomial(i -> 3 - float(i), monomial_vector([x * x, x]))
    for p in [
        polynomial(i -> float(i), [x, x * x]),
        polynomial(i -> 1.0, [x * x, x, x * x]),
        polynomial(i -> float(i), monomial_vector([x * x, x])),
    ]
        @test collect(coefficients(p)) == [1.0, 2.0]
        @test collect(monomials(p)) == monomial_vector([x, x^2])
    end

    @test (x + y)' == x + y
    @test transpose(x + y) == x + y
    @test transpose([1 2; 3 4] * x) == [1 3; 2 4] * x

    p = 4x^2 * y + x * y + 2x
    @test remove_monomials(p, [x * y]) == 4x^2 * y + 2x
    @test filter_terms(OfDegree(2), 4x^2 * y + x * y + 2x) == x * y
    @test filter_terms(!OfDegree(2), 4x^2 * y + x * y + 2x) == 4x^2 * y + 2x
    @test filter_terms(iseven ∘ MP.coefficient, 4x^2 * y + x * y + 2x) ==
          4x^2 * y + 2x

    @test_throws MethodError push!([1], x + 1)

    @test polynomial([1 2; 3 4], [x^2, y]) == x^4 + 5x^2 * y + 4y^2
    @test polynomial([1 2; 3 4], [x^2, y], Float64) isa
          AbstractPolynomial{Float64}
    @test polynomial([1 2; 3 4], [y, x^2]) == y^2 + 5x^2 * y + 4x^4
    @test polynomial([1 2; 3 4], [y, x^2], Float64) isa
          AbstractPolynomial{Float64}
    @test polynomial([1 2; 3 4], monomial_vector([y, x^2])) ==
          4x^4 + 5x^2 * y + y^2
    @test polynomial([1 2; 3 4], monomial_vector([y, x^2]), Float64) isa
          AbstractPolynomial{Float64}

    @test (@inferred round(2.6x + 1.001x^2)) == 3x + 1x^2
    @test (@inferred round(3.1x * y)) == 3x * y
    @test (@inferred round(2.613x + 1.1051x^2, digits = 2)) ≈ 2.61x + 1.11x^2
    @test (@inferred round(3.145x * y, digits = 1)) ≈ 3.1x * y

    @testset "Graded Lex Order" begin
        Mod.@polyvar x y z
        p = 3 * y^2 + 2 * y * x
        @test collect(coefficients(p)) == [3, 2]
        @test collect(monomials(p)) == monomial_vector([x * y, y^2])
        # Examples from p. 59 of the 4th edition of "Ideals, Varieties, and Algorithms" of Cox, Little and O'Shea
        f = 4 * x * y^2 * z + 4 * z^2 - 5 * x^3 + 7 * x^2 * z^2
        @test collect(coefficients(f)) == [4, -5, 4, 7]
        @test collect(monomials(f)) ==
              monomial_vector([x^2 * z^2, x * y^2 * z, x^3, z^2])
        #@test ordering(f) === GradedLex #FIXME
    end

    @testset "Convertion" begin
        Mod.@polyvar x y z
        p = 2.5x + 1 - 2.5x
        @test MP.convert_to_constant(Int, p) == 1
        @test convert(typeof(p), p) === p
        @test convert(Union{Nothing,typeof(p)}, p) === p
        a = 2y
        q = polynomial([a, z, -a], [x, 1, x])
        r = @inferred convert_to_constant(q)
        @test r == z
        @test MA.add!!(r, 1) == z + 1
        @test convert_to_constant(q) == z
        @test iszero(@inferred convert_to_constant(zero(q)))
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

    @testset "$f is a mutable copy, see issue DynamicPolynomials#62" for f in [
        zero,
        one,
    ]
        p = 2x + 1 + zero(x + y)
        q = f(p)
        q = MA.add!!(q, 2y)
        @test q == f(p) + 2y
        @test p == 2x + 1
        q = f(2x + 1)
        mutable =
            MA.mutability(typeof(q), +, typeof(q), typeof(2y)) isa MA.IsMutable
        result = @inferred MA.add!!(q, 2y)
        @test result == f(2x + 1) + 2y
        @test (result === q) == mutable
        @test q == (mutable ? result : f(2x + 1))
    end

    @testset "map_coefficients $nz" for nz in [false, true]
        p = 2x + 1
        @test map_coefficients(x -> x / 2, p, nonzero = nz) == 1.0x + 0.5
        @test map_coefficients(x -> x / 2, 3x, nonzero = nz) == 1.5x
        @test p === SA.map_coefficients!(x -> x + 1, p, nonzero = nz)
        @test p == 3x + 2
        q = zero(p)
        @test q === SA.map_coefficients_to!(q, x -> 2x, p, nonzero = nz)
        @test q == 6x + 4
        @test q === SA.map_coefficients_to!(q, x -> 2x, 3x, nonzero = nz)
        @test q == 6x
        @test q === SA.map_coefficients_to!(q, x -> 2x, x, nonzero = nz)
        @test q == 2x
        @test q === SA.map_coefficients_to!(q, x -> 2x, x^2, nonzero = nz)
        @test q == 2x^2

        input = 3x^2
        output = zero(x + 1)
        SA.map_coefficients_to!(output, identity, input; nonzero = nz)
        @test output == 3x^2
        if exponents(leading_monomial(output)) isa AbstractVector
            exponents(leading_monomial(output))[1] = 3
            @test output == 3x^3
        end
        @test input == 3x^2

        for input in (x, x^2, 3x, 2x + 1)
            output = one(x + y)
            @test (@inferred SA.map_coefficients_to!(
                output,
                c -> 2c,
                input;
                nonzero = nz,
            )) === output
            @test output == 2input
            @test variables(output) == variables(x + y)
        end
        for input in (y, y^2, 3y, 2y + 1)
            output = x + 1
            if MA.promote_operation(+, typeof(output), typeof(input)) ===
               typeof(output)
                @test (@inferred SA.map_coefficients_to!(
                    output,
                    identity,
                    input;
                    nonzero = nz,
                )) === output
                @test output == input
                @test variables(output) == variables(x + y)
            else
                @test_throws ArgumentError SA.map_coefficients_to!(
                    output,
                    zero,
                    input;
                    nonzero = nz,
                )
                @test output == x + 1
            end
        end
    end
end

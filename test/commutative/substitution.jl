if !@isdefined(Mod)
    include("../SimplePolynomials.jl")
    import .SimplePolynomials as Mod
    include("../utils.jl")
end

using Test
import Test: @inferred

import MutableArithmetics as MA

function substitution_allocated(p, vals)
    p(vals)
    return @allocated p(vals)
end

@testset "Substitution" begin
    Mod.@polyvar x[1:3]

    @test subs(2, x[1] => 3) == 2
    @test iszero((x[1] - x[1])(x[1] => x[2]))

    a = (x[1])(x[1] => x[2])
    b = x[2]
    @test (x[1])(x[1] => x[2]) == x[2]

    p = x[1] * x[2] * x[3]
    @test p(x => (1, 2, 3)) == 6

    p = x[1]^2 + x[1] * x[3] - 3
    @test MP.convert_to_constant(Int, p(x => (5, x[1] * x[2], 4))) == 42

    p = x[1]^2 + x[2]^2
    q = p(x[1:2] => [1 -1; 1 1] * vec(x[1:2]))
    @test q == 2p

    t = 2.0 * x[1] * x[2]^2
    @test t(x => [1, 2, 3]) == 8.0
    @test t((x[2], x[1]) => [1, 3]) == 6.0

    p = x[1] + x[2] + 2 * x[1]^2 + 3 * x[1] * x[2]^2
    @test (@inferred p((x[1], x[2]) => (1.0, 2.0))) == 17.0

    m = x[1] * x[2]
    @inferred subs(m, x[2] => 2.0)
    @inferred subs(p, x[2] => 2.0)
    @test subs(p, x[2] => 2.0) == 13x[1] + 2 + 2x[1]^2
    @inferred subs(p, x[2] => x[1])
    @inferred subs(p, x[2] => x[1]) == 2x[1] + 2x[1]^2 + 3x[1]^3

    p = x[1] - x[1]
    @test iszero(p)
    @inferred p(x[1] => 1)
    @test iszero(p(x[1] => 1))
    @inferred subs(p, x[1] => 1)
    @test iszero(subs(p, x[1] => 1))

    # https://github.com/JuliaAlgebra/DynamicPolynomials.jl/issues/141
    @test subs(x[1]^2 + x[2]^2, x[1] => x[3]^2) == x[2]^2 + x[3]^4

    q = (x[1] + 1) / (x[1] + 2)
    @test isapproxzero(q(x[1] => -1))
    @test !isapproxzero(q(x[1] => 1))
    @test q(x[1] => 1) ≈ 2 / 3

    q = (x[1] + x[3]) / (x[2] - 1)
    @test subs(q, x[1:2] => (x[2], x[1])) == (x[2] + x[3]) / (x[1] - 1)

    # Taken from perf/runbenchmark.jl
    p = x[1] + x[2] + 2x[1]^2 + x[2]^3

    varst = (x[1], x[2])
    varsv = [x[1], x[2]]
    valst = (1.0, 2.0)
    valsv = [1.0, 2.0]

    @test variable.(Tuple(varsv)) == varst
    @test variable.(Tuple(varsv)) isa
          Tuple{<:AbstractVariable,<:AbstractVariable}

    @test p(variables(p) => valst) == 13
    @test p(variables(p) => valsv) == 13
    @test p(varst => valst) == 13
    @test p(varsv => valsv) == 13

    @testset "Positional evaluation" begin
        for q in (x[1], x[1] * x[2], 2x[1] * x[2], p)
            vals = ntuple(i -> Float64(i), nvariables(q))
            expected = q(variables(q) => vals)
            @test (@inferred q(vals)) == expected
            @test (@inferred q(vals...)) == expected
            @test (@inferred q(collect(vals))) == expected
            @test (@inferred substitute(Eval(), q, collect(vals))) == expected
            @test (@inferred substitute(Eval(), q, vals)) == expected
            @test_throws ArgumentError q([vals...; 3.0])
            @test_throws ArgumentError q(Float64[])
        end
        c = constant_monomial(typeof(x[1]))
        @test (@inferred c(ntuple(_ -> 2, nvariables(c)))) == 1
    end

    @testset "Mutable substitution values" begin
        # https://github.com/JuliaAlgebra/DynamicPolynomials.jl/issues/208
        Mod.@polyvar a b
        m = a * b
        p = a * b * (1 + b)
        for T in (BigInt, BigFloat, Rational{BigInt})
            # The two entries deliberately share their mutable storage.
            vals = fill(T(-1), 2)
            original = MA.mutable_copy(vals)
            for point in (vals, Tuple(vals))
                @test (@inferred m(point)) == 1
                @test vals == original
                @test iszero(@inferred p(point))
                @test vals == original
            end
            @test iszero(@inferred p(a => vals[1], b => vals[2]))
            @test vals == original
            @test iszero(@inferred subs(p, a => vals[1], b => vals[2]))
            @test vals == original
            @test (@inferred subs(p, a => vals[1])) == -b - b^2
            @test vals == original
        end
    end

    @testset "Mutable source coefficients" begin
        # https://github.com/JuliaAlgebra/DynamicPolynomials.jl/pull/100
        Mod.@polyvar a b
        for T in (BigInt, BigFloat, Rational{BigInt})
            t = term(T(2), a)
            q = @inferred subs(t, a => b)
            @test q == 2b
            MA.operate!(one, coefficient(q))
            @test coefficient(t) == 2

            for value in (false, true)
                t = term(T(2), a)
                result = @inferred t(a => value)
                @test result == 2value
                MA.operate!(one, result)
                @test coefficient(t) == 2

                p = T(2) * a + T(3) * b
                q = @inferred subs(p, a => value)
                @test q == 2value + 3b
                MA.operate!(one, coefficient(q, b))
                @test p == 2a + 3b

                # Combining equal remaining monomials must also preserve p.
                p = T(2) * a * b + T(3) * b
                q = @inferred subs(p, a => value)
                @test q == (2value + 3) * b
                @test p == 2a * b + 3b
                MA.operate!(one, coefficient(q, b))
                @test p == 2a * b + 3b
            end

            p = T(2) * a + T(3)
            q = @inferred subs(p, a => b)
            @test q == 2b + 3
            MA.operate!(one, coefficient(q, b))
            @test p == 2a + 3
        end
    end

    @testset "Heterogeneous tuple substitution" begin
        Mod.@polyvar a b c d e f g h
        m = a * b * c * d * e * f * g * h
        @test (@inferred subs(m, a => 2.0)) == 2.0b * c * d * e * f * g * h
        @test subs(m + 1, a => 2.0) == 2.0b * c * d * e * f * g * h + 1
        for vals in (fill(2, 8), ntuple(_ -> 2, 8))
            @test substitution_allocated(m, vals) == 0
        end
        # Powers allocate BigInts, but multiplication reuses the accumulator.
        # This budget excludes allocating another product for every factor.
        @test substitution_allocated(m, fill(big(2), 8)) <= 16 * BIGINT_ALLOC
    end

    @testset "Eval on constant monomial" begin
        m = constant_monomial(prod(x))
        mr = m(x => reverse(x))
        @test isconstant(mr)
        @test variables(mr) == x
    end

    @testset "Partial numeric substitution" begin
        Mod.@polyvar a b c
        q = 3a^2 * b + 2b * c + c^3 + 1
        r = @inferred subs(q, a => 2)
        @test r == 12b + 2b * c + c^3 + 1
        @test variables(r) == variables(b * c)
        @test (@inferred subs(q, a => 2.0, c => 1)) == 14.0b + 2.0
        @test (@inferred subs(q, (a, b) => (1, 0.5))) == c^3 + c + 2.5
        @test (@inferred subs(q, b => 1 + im)) ==
              (3 + 3im) * a^2 + (2 + 2im) * c + c^3 + 1
        # Terms with the same remaining monomial are combined
        @test (@inferred subs(a * b - 2b + c, a => 2)) == c
    end

    Mod.@polyvar x y
    @test subs([x^2 + y, x + y], x => y) == [y^2 + y, 2y]

    @testset "Typed replacement values" begin
        Mod.@polyvar z
        m = x * y
        @test (@inferred MA.promote_operation(^, typeof(2x), Int)) ===
              typeof(2x)
        @test (@inferred subs(m, y => 2.0)) isa AbstractTerm{Float64}
        @test (@inferred subs(m, y => 2.0)) == 2.0x
        @test (@inferred subs(m, x => 3, y => 2.0)) == 6.0
        @test (@inferred m(x => 3, y => 2.0)) === 6.0
        @test (@inferred m([x, y] => (3, 2.0))) === 6.0
        @test (@inferred m((x, y) => [3.0, 2.0])) === 6.0
        @test (@inferred m([x, y] => [3.0, 2.0])) === 6.0
        @test (@inferred m([x^1, y^1] => [3.0, 2.0])) === 6.0
        @test (@inferred subs(m, (x,) => (z + 1,), y => 2.0)) == 2z + 2
        @test (@inferred subs(m, x => z, (y,) => (2.0,))) == 2.0z
        @test (@inferred subs(x + 2y, x => y, y => x)) == y + 2x
        @test (@inferred subs(m)) == m
        @test (@inferred substitute(Eval(), m, (x, y) => (3, 2.0))) === 6.0
        @test (@inferred substitute(Eval(), m, x => 3, y => 2.0)) === 6.0
        @test_throws ArgumentError m(x => 2)
        @test_throws ArgumentError subs(m, (x, y) => [1])
        @test_throws ArgumentError subs(m, [x] => [1, 2])
        @test_throws ArgumentError subs(m, [x, y] => (1,))

        for p in (x, m, 2m, m + x + 1),
            s in ((x => 2.0, y => 3), ((x, y) => (z, z + 1),)),
            st in (Subs(), Eval())

            q = @inferred substitute(st, p, s)
            @test (@inferred MA.promote_operation(
                substitute,
                typeof(st),
                typeof(p),
                typeof.(s)...,
            )) === typeof(q)
        end

        p = x^2 + 2y + 1
        @test (@inferred subs(p, x => z + 1.0, y => 2)) == z^2 + 2z + 6
        @test p == x^2 + 2y + 1
    end

    @testset "Subs with no variables" begin
        Mod.@polyvar x
        t = convert(term_type(x, Int), 3)
        @test t == @inferred subs(t, x => x + 1)
        @test t == @inferred subs(t, x => x + 1.0)
        @test t == @inferred subs(t, x => 1x)
        @test t == @inferred subs(t, x => 1.0x)
        @test t == @inferred subs(t, x => 1.0)
        @test (@inferred t()) === 3
        @test (@inferred subs(t)) == t
        m = constant_monomial(typeof(x))
        @test (@inferred m()) === 1
        @test isone(@inferred subs(m))
        @test iszero(@inferred zero(t)(x => 2.0))
    end

    for (p, s) in [(x^1, x => 2x)]
        @test MA.promote_operation(substitute, Subs, typeof(p), typeof(s)) ==
              typeof(subs(p, s))
    end

    @testset "Issue #302 $T" for T in [BigFloat, BigInt]
        Mod.@polyvar x[1:3]
        t = T(1) * x[1]
        @test copy(t).coefficient !== t.coefficient
        @test MA.mutable_copy(t).coefficient !== t.coefficient
        F = T(5) * x[1] * x[2] * x[3] + T(1) * x[1] * x[2]^2
        @test subs(F, x[3] => T(0)) == x[1] * x[2]^2
        @test subs(F, x[3] => 0) == x[1] * x[2]^2
    end
end

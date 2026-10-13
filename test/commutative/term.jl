if !@isdefined(Mod)
    import DynamicPolynomials as Mod
    include("../utils.jl")
end

import MutableArithmetics as MA
import MultivariatePolynomials as MP

struct CoefNotComparable end
Base.iszero(::CoefNotComparable) = false

@testset "Term" begin
    Mod.@polyvar x
    @test convert(Any, 1x) == 1x
    @test one(1x) == one(1.0x) == 1
    @test zero(1x) == zero(1.0x) == 0
    @test nvariables(0.0x) == 1
    @test nvariables(1x) == 1
    #@inferred one(1x)
    @inferred zero(1x)
    #@inferred one(1.0x)
    @inferred zero(1.0x)

    @test !isapproxzero((1 + 1e-10im) * x^2)
    @test isapproxzero(1e-10im * x^2)

    @test monic(2.0x) isa AbstractTerm{Float64}
    @test monic(2.0x) == x
    @test monic(2x^2) isa AbstractTerm{Int}
    @test monic(2x^2) == x^2

    @test leading_term(2x^2) == 2x^2
    @test nterms(2x^2) == 1
    @test terms(2x^2) == [2x^2]
    @test nterms(0 * x) == 0
    @test terms(0 * x) == typeof(0 * x)[]
    @test nterms(0.0x) == 0
    @test terms(0.0x) == typeof(0.0x)[]

    @test nterms(polynomial(0.0x)) == 0
    @test nterms(convert(polynomial_type(0.0x), 0.0x)) == 0

    @test term(x) isa AbstractTerm
    @test term(x^2) == x^2
    @test term(1x^2) isa AbstractTerm
    @test term(1x) == x
    @test zero_term(1x) == 0 * x

    Mod.@polyvar y
    @testset "Products of terms stay terms" begin
        for (a, b, expected) in (
            (2x, 3y, term(6, x * y)),
            (2x, 3.0y, term(6.0, x * y)),
            (big(2) * x, big(3) * y, term(big(6), x * y)),
            (0 * x, 3y, term(0, x * y)),
        )
            product = @inferred a * b
            @test product isa typeof(expected)
            @test product == expected
            @test (@inferred MA.promote_operation(*, typeof(a), typeof(b))) ===
                  typeof(product)
            @test monomial(a) == x
            @test monomial(b) == y
            if coefficient(product) isa BigInt
                MA.operate!(+, coefficient(product), big(1))
                @test coefficient(a) == 2
                @test coefficient(b) == 3
            end
        end
    end
    @testset "Monomial products preserve terms" begin
        for c in (0, 2, 2.5, 2 // 3, big(2))
            t = term(c, y)
            for m in (x, x^2), (a, b) in ((m, t), (t, m))
                product = @inferred a * b
                @test product isa AbstractTerm{typeof(c)}
                @test product == term(c, m * y)
                @test (@inferred MA.promote_operation(
                    *,
                    typeof(a),
                    typeof(b),
                )) === typeof(product)
                @test (@inferred MA.operate!!(*, a, b)) == product
                @test monomial(t) == y
                if c isa BigInt
                    MA.operate!(+, coefficient(product), big(1))
                    @test coefficient(t) == 2
                end
            end
        end
    end
    @test degree(2x^2, x) == 2
    @test degree(2x^2, y) == 0
    @test degree(2x^2, y) == 0

    t = 3x^2 * y^4
    alloc_test(() -> convert(typeof(t), t), 0)
    typetests(t)
    typetests([t, 2x])
    @test (@inferred polynomial(t)) isa AbstractPolynomial{Int}
    @test (@inferred polynomial(t, Float64)) isa AbstractPolynomial{Float64}

    @test_throws InexactError push!([1], 2x)
    @test_throws InexactError push!([x^2], 2x)

    t = MP.term(1, x)
    @test convert(term_type(t, Number), t) isa term_type(t, Number)

    @testset "Effective variables" begin
        T = variable_union_type(x)
        @test x isa T
        @test y isa T
        @test T[x, y] == @inferred effective_variables(x * y)
        @test T[x, y] == @inferred effective_variables(y * x)
        @test T[x] == @inferred effective_variables(x * y^0)
        @test T[x] == @inferred effective_variables(y^0 * x)
        @test T[y] == @inferred effective_variables(x^0 * y)
        @test T[y] == @inferred effective_variables(y * x^0)
    end

    @testset "Compare terms" begin
        t1 = 1 * x
        t2 = 2 * x
        @test t1 < t2
        @test t1 <= t2
        @test !(t1 > t2)
        @test !(t1 >= t2)

        t1 = (1 + 1im) * x
        t2 = (2 + 2im) * x
        @test t1 < t2
        @test t1 <= t2
        @test !(t1 > t2)
        @test !(t1 >= t2)

        a = CoefNotComparable()
        b = CoefNotComparable()
        t1 = term(a, x)
        t2 = term(b, y)
        @test t1 > t2
        @test t1 >= t2
        @test !(t1 < t2)
        @test !(t1 <= t2)
        t1 = term(a, x)
        t2 = term(b, x)
        @test !(t1 > t2)
        @test t1 >= t2
        @test !(t1 < t2)
        @test t1 <= t2
    end

    @testset "MA $T" for T in [Int, BigInt]
        t = one(T) * x
        s = @inferred MA.operate!!(*, t, x)
        @test monomial(s) == x^2
        @test t == one(T) * x
        u = @inferred MA.operate!!(*, s, term(T(2), x^3))
        @test u isa AbstractTerm{T}
        @test u == term(T(2), x^5)
        @test s == term(one(T), x^2)
        @test t == one(T) * x
    end
end

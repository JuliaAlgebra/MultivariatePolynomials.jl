import MutableArithmetics as MA
import MultivariatePolynomials as MP

struct CoefNotComparable end
Base.iszero(::CoefNotComparable) = false

struct Term2{T,M} <: MP.AbstractTermLike{T}
    monomial::M
end
MP.coefficient(t::Term2{T}) where {T} = 2one(T)
MP.monomial(t) = t.monomial
MP.term_type(::Type{Term2{T,M}}) where {T,M} = MP.Term{T,M}

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
    @test convert(MP.Term{Number,MP.monomial_type(t)}, t) isa
          MP.Term{Number,MP.monomial_type(t)}

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
        t1 = a * x
        t2 = b * y
        @test t1 > t2
        @test t1 >= t2
        @test !(t1 < t2)
        @test !(t1 <= t2)
        t1 = a * x
        t2 = b * x
        @test !(t1 > t2)
        @test t1 >= t2
        @test !(t1 < t2)
        @test t1 <= t2
    end

    @testset "MA $T" for T in [Int, BigInt]
        M = typeof(x^2)
        t = one(T) * x
        s = MA.operate!!(*, t, x)
        @test monomial(s) == x^2
        if T == BigInt && MA.mutability(M) isa MA.IsMutable
            @test monomial(t) == x^2
        end
        u = MA.operate!!(*, s, Term2{T,M}(x^3))
        @test monomial(u) == x^5
        @test coefficient(u) == 2
        if T == BigInt && MA.mutability(M) isa MA.IsMutable
            @test monomial(t) == x^5
            @test coefficient(t) == 2
        end
    end
end

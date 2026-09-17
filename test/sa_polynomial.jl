module TestSAPolynomial

using Test
import DynamicPolynomials as DP
import MultivariatePolynomials as MP
import StarAlgebras as SA
import LinearAlgebra

@testset "Shared polynomial representation" begin
    DP.@polyvar x y
    p = 2x + y
    @test x isa MP.AbstractPolynomialLike
    @test x^2 isa MP.AbstractTermLike
    @test 2x isa MP.AbstractTerm
    @test p isa MP.AbstractPolynomial
    @test typeof(p) === DP.Polynomial{Int}
    @test MP.term_type(x) === DP.Term{Int}
    @test MP.polynomial_type(x) === typeof(p)
    @test MP.polynomial_type(p, Float64) === DP.Polynomial{Float64}
    @test MP.coefficient_type(x^2) === Int
    @test MP.coefficients(p) == [1, 2]
    @test MP.coefficient.(MP.terms(p)) == [1, 2]
    @test length(MP.monomials(p)) == 2
    @test MP.polynomial(p) === p
    @test MP.polynomial!(p) === p
    @test MP.polynomial!(2x, Float64) isa DP.Polynomial{Float64}
    @test sprint(show, typeof(p)) == "DynamicPolynomials.Polynomial{$Int}"
    @test sprint(show, typeof(2x)) == "DynamicPolynomials.Term{$Int}"
    @test promote_type(typeof(x), String) === Any
    @test_throws MethodError p * "x"
    @test_throws InexactError p * 0.5
    @test MP.coefficients(p / 2) == [0.5, 1.0]
    @test MP.coefficients(p // 2) == [1 // 2, 1 // 1]
    @test div(p, 2) == MP.polynomial(x)
    @test_throws InexactError p / 2.5
    @test iszero(p - (x + y) - x)
    @test iszero((x + y)^2 - (x^2 + 2x*y + y^2))
    @test iszero(sum([x, y]) - (x + y))
    @test iszero(sum([p, p]) - 2p)
    @test MP.coefficients(p) == [1, 2]
    @test iszero(sum(DP.Polynomial{Int}[]))
    @test !isone(zero(p))
    @test SA.coeffs(zero(p)) != SA.coeffs(p)
    @test MP.coefficients(MP.differentiate(p, x)) == [2]

    @test MP.coefficient(p, x) == 2
    @test MP.coefficient(x, y) == 0
    @test MP.coefficients(x, [x, y]) == [1, 0]
    @test MP.MA.scaling(one(x)) == 1
    @test p - y == MP.polynomial(2x)
    @test MP.subs(p, x => 2, y => 3) == 7
    @test iszero(MP.subs(p, x => 2) - (y + 4))
    @test MP.monomial([x, y], [1, 0]) == MP.monomial(x)
    @test hash(MP.monomial([x, y], [1, 0])) == hash(MP.monomial(x))

    @test p(x => 2, y => 3) == 7
    @test (2x)(x => 3) == 6
    @test MP.coefficient_type([x, y]) === Int

    DP.@polyvar u v monomial_order=MP.LexOrder
    q = u + v
    @test !(q isa DP.Polynomial)
    @test MP.polynomial_type(u) === typeof(q)
    @test iszero(sum([u, v]) - q)
end

@testset "Equality and adjoints" begin
    DP.@polyvar x y
    p = x + y
    for comp in (==, isequal)
        @test comp(zero(p), 0) && comp(0, zero(p))
        @test comp(one(p), 1) && comp(1, one(p))
        @test comp(MP.polynomial(x), x) && comp(x, MP.polynomial(x))
        @test comp(MP.term(2, one(x)), 2) && comp(2, MP.term(2, one(x)))
        @test comp(one(x), 1) && comp(1, one(x))
        @test !comp(p, 0) && !comp(0, p)
    end
    @test adjoint(p) == transpose(p) == p
    constant = MP.term(2, one(x))
    @test hash(constant) == hash(MP.polynomial(constant)) == hash(2)
    @test haskey(Dict{Any,Bool}(constant => true), 2)
    t = (1 + 2im) * x
    @test adjoint(t) == (1 - 2im) * x
    @test transpose(t) == t
    c = t + (3 - im) * y
    @test adjoint(c) == (1 - 2im) * x + (3 + im) * y
    @test transpose(c) == c
    @test LinearAlgebra.dot(p, p) == p * p
    @test LinearAlgebra.dot(t, t) == 5 * x^2
    @test LinearAlgebra.dot([p, p], [p, p]) == 2 * p * p
    @test LinearAlgebra.dot([t, t], [t, t]) == 10 * x^2
end

@testset "Matching coefficient type" begin
    DP.@polyvar x
    a = [1 2; 3 4]
    t = MP.term(a, x)
    c = one(x)
    for p in (t + a, a + t)
        @test MP.coefficient(p, x) == a
        @test MP.coefficient(p, c) == a
    end
    @test MP.coefficient(t - a, c) == -a
    @test MP.coefficient(a - t, x) == -a
    @test MP.term(a, c) == a
    @test MP.polynomial(MP.term(a, c)) == a
end

@testset "Coefficient types in division and GCD" begin
    DP.@polyvar x y
    algo = MP.SubresultantAlgorithm()
    immutable = MP.MA.IsNotMutable()
    for (p, T, field) in (
        (x, Int, Rational{Int}),
        (x^2, Int, Rational{Int}),
        (2x, Int, Rational{Int}),
        (2x + 1, Int, Rational{Int}),
        ((2 // 3) * x, Rational{Int}, Rational{Int}),
        ((2 // 3) * x + 1, Rational{Int}, Rational{Int}),
        (2.0x, Float64, Float64),
        (2.0x + 1, Float64, Float64),
    )
        P = typeof(p)
        for op in (div, rem)
            @test MP.MA.promote_operation(op, P, P) === DP.Polynomial{field}
        end
        for op in (MP.pseudo_rem, MP.rem_or_pseudo_rem)
            @test MP.MA.promote_operation(op, P, P, typeof(algo)) ===
                  DP.Polynomial{T}
        end
    end
    for p in (x, x^2)
        @test (@inferred MP.content(p, algo, immutable)) === 1
        @test MP.primitive_part(p, algo, immutable) == p
    end
    @test MP.content(6x, algo, immutable) == 6
    @test MP.content(6x + 9, algo, immutable) == 3
    @test MP.content(zero(x + y), algo, immutable) == 0
    @test MP.primitive_part(6x + 9, algo, immutable) == 2x + 3
    for p in (2.5x, 2.5x + 1)
        @test MP.content(p, algo, immutable) === 1.0
        @test MP.primitive_part(p, algo, immutable) === p
    end
    complex_term = (2.0 + 3.0im) * x
    @test MP.primitive_part(complex_term, algo, immutable) === complex_term
    @test gcd(x^2, x^3) == x^2
    @test gcd(6x, 9x) == 3x
    @test isempty(Test.detect_unbound_args(MP))
end

@testset "Monomial divisibility" begin
    DP.@polyvar x y z
    for (a, b, expected) in (
        (x, x, true),
        (x, y, false),
        (x, x^2, true),
        (x^2, x, false),
        (x^2, x^2, true),
        (x^2, y^2, false),
        (x * y, x^2 * y^3, true),
        (x * y^3, x^2 * y, false),
        (y, x * y * z, true),
        (x * z, y^2, false),
        (x * y, x^2, false),
        (y * z, y^2, false),
        (one(x), x^2, true),
        (x^2, one(x), false),
        (one(x), one(y), true),
        (MP.monomial([x, y, z], [0, 2, 0]), y^3, true),
        (y, MP.monomial([x, y, z], [0, 2, 0]), true),
        (MP.term(3, x * y), MP.term(2, x^2 * y), true),
        (MP.term(2, x^2 * y), MP.term(3, x * y), false),
    )
        originals = deepcopy((a, b))
        @test MP.divides(a, b) === expected
        @test (a, b) == originals
    end
    @test (@inferred MP.divides(x, x^2))
    DP.@ncpolyvar u v
    @test MP.divides(u, u)
    @test !MP.divides(u, v)
    @test_throws ErrorException MP.divides(u, u^2)
    @test_throws ErrorException MP.divides(u^2, v^2)
end

@testset "Polynomial conversion" begin
    DP.@polyvar x
    R = DP.Polynomial{Rational{Int}}
    for p in (x, x^2, 2x, 2x + 1, 3, 0)
        converted = convert(R, p)
        @test converted isa R
        @test converted == p
        @test convert(R, converted) === converted
    end
    p = 2x + 1
    @test parent(convert(R, p)) === parent(p)
    @test_throws InexactError convert(DP.Polynomial{Int}, (1 // 2) * x + 1)
    @test_throws InexactError convert(DP.Polynomial{Int}, (1 // 2) * x)
    @test_throws InexactError convert(DP.Polynomial{Int}, 1 // 2)
end

@testset "Polynomial conversion in division and GCD" begin
    DP.@polyvar x y
    for p in (x, x^2, 2x, 2x + 1)
        original = deepcopy(p)
        q, r = divrem(p, typeof(p)[])
        @test isempty(q)
        @test r == p
        @test r isa DP.Polynomial{Rational{Int}}
        MP.MA.operate!(zero, r)
        @test p == original
    end
    q, r = divrem(zero(2x), 2x)
    @test q isa DP.Polynomial{Rational{Int}}
    @test iszero(q) && iszero(r)

    algo = MP.GeneralizedEuclideanAlgorithm()
    g = MP.primitive_univariate_gcd!(2x, zero(2x), algo)
    @test g isa DP.Polynomial{Int}
    @test g == 2x
    @test MP.inflate(3, one(x), one(x)) == 3
    nested = MP.term(2x, MP.monomial(y))
    c = MP.content(nested, algo, MP.MA.IsNotMutable())
    @test c isa DP.Term{Int}
    @test c == 2x
end

@testset "GCD result type promotion" begin
    DP.@polyvar x y
    for (a, b) in (
        (x, y),
        (x^2, x),
        (6x^2, 9x),
        (2x, y),
        (2.0x, 3y),
        ((2 // 3) * x, 3y),
        ((2 + 3im) * x, 3y),
        ((2 + 3im) * x, 3.0y),
        (MP.term(2x, MP.monomial(y)), MP.term(3x, MP.monomial(y))),
    )
        for (p, q) in ((a, b), (b, a))
            @test MP.MA.promote_operation(gcd, typeof(p), typeof(q)) ===
                  typeof(gcd(p, q))
            for algo in
                (MP.GeneralizedEuclideanAlgorithm(), MP.SubresultantAlgorithm())
                @test MP.MA.promote_operation(
                    gcd,
                    typeof(p),
                    typeof(q),
                    typeof(algo),
                ) === typeof(gcd(p, q, algo))
            end
        end
    end
    @test MP.MA.promote_operation(gcd, typeof(x + 1), typeof(x + 1)) ===
          DP.Polynomial{Int}
    nested = MP.term(x, MP.monomial(y))
    c = MP.content(nested, MP.SubresultantAlgorithm(), MP.MA.IsNotMutable())
    @test c isa MP.AbstractMonomial
    @test c == x
end

end # module

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

end # module

if !@isdefined(Mod)
    import DynamicPolynomials as Mod
    include("../utils.jl")
end

import MutableArithmetics as MA
import StarAlgebras as SA

@testset "Noncommutative term products" begin
    Mod.@ncpolyvar x y
    @testset "Monomial products preserve terms" begin
        for (a, b, expected) in (
            (x, 2y, term(2, x * y)),
            (2y, x, term(2, y * x)),
            (x * y, 3x^2, term(3, x * y * x^2)),
            (3x^2, x * y, term(3, x^3 * y)),
        )
            product = @inferred a * b
            @test product isa AbstractTerm{Int}
            @test product == expected
            @test (@inferred MA.promote_operation(*, typeof(a), typeof(b))) ===
                  typeof(product)
            @test (@inferred MA.operate!!(*, a, b)) == expected
        end
    end
    @testset "Products of terms are polynomials" begin
        for (a, b, expected) in (
            (2x, 3y, polynomial(term(6, x * y))),
            (2x, 3x, polynomial(term(6, x^2))),
        )
            product = @inferred a * b
            @test product isa AbstractPolynomial{Int}
            @test product == expected
            @test (@inferred MA.promote_operation(*, typeof(a), typeof(b))) ===
                  typeof(product)
        end
        @test (@inferred (3y) * (2x)) == polynomial(term(6, y * x))
    end
end

@testset "Noncommutative monomial keys" begin
    Mod.@ncpolyvar x y z
    b = @inferred MP.FullBasis{MP.Monomial}([x, y, x])
    @test variables(b) == [x, y]
    @test nvariables(b) == 2
    @test SA.key_type(b) === monomial_type(x)
    @test !(SA.object(MP.algebra(b)) isa MP.Variables)
    @test (@inferred b[x*y*x]) == x * y * x
    @test haskey(b, x * y * x)
    @test !haskey(b, x * z)
    @test_throws ArgumentError b[x*z]
    @test (@inferred collect(Iterators.take(b, 7))) == monomials([x, y], 0:2)
    empty_basis = MP.FullBasis{MP.Monomial}(typeof(x)[])
    @test collect(empty_basis) == [one(x)]
    @test !isempty(empty_basis)
    @test b == MP.FullBasis{MP.Monomial}([y, x])
    @test hash(b) == hash(MP.FullBasis{MP.Monomial}([y, x]))
    @test_throws ErrorException MP.FullBasis{MP.AbstractMonomialBasis}([x, y])

    p = @inferred 2x * y * x + 3y * x * y
    @test SA.key_type(SA.coeffs(p)) === monomial_type(x)
    @test SA.keys(SA.coeffs(p)) == monomials(p)
    @test (@inferred p * (x + y)) ==
          2x * y * x^2 + 2x * y * x * y + 3y * x * y * x + 3y * x * y^2
    @test (@inferred (x + y) * p) ==
          2x^2 * y * x + 3x * y * x * y + 2y * x * y * x + 3y^2 * x * y
    @test (p * (x + 1)) * (y + 2) == p * ((x + 1) * (y + 2))
    @test (@inferred p(x => 2, y => 3)) == 78
    @test (@inferred subs(p, x => y * x)) == 2y * x * y^2 * x + 3y^2 * x * y
    @test zero(p)(x => 2, y => 3) == 0
    @test MP.implicit(p) == p
    @test variables(zero(p)) == [x, y]
    @test one(p) == 1
    @test iszero(zero(typeof(p)))
    @test iszero(zero(x * y))
    @test iszero(zero(typeof(x * y)))
    @test (@inferred MP.polynomial_type(b, Float64)) ===
          polynomial_type(x, Float64)

    q = z + 1
    pp, qq = @inferred SA.promote_bases(p, q)
    @test parent(pp) == parent(qq)
    @test variables(pp) == [x, y, z]
    @test SA.coeffs(pp) === SA.coeffs(p)
    @test pp == p
    @test qq == q
    output = zero(p)
    @test (@inferred MA.operate!(+, output, q)) === output
    @test output == q
    @test variables(output) == [x, y, z]
    @test variables(p) == [x, y]
end

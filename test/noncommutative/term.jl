if !@isdefined(Mod)
    import DynamicPolynomials as Mod
    include("../utils.jl")
end

import MutableArithmetics as MA

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
        # Algebra elements currently have a fixed variable sequence: the common
        # basis for x and y cannot represent the reversed word y*x.
        @test_broken (3y) * (2x) == polynomial(term(6, y * x))
    end
end

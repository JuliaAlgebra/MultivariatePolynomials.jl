using LinearAlgebra, Test
using Combinatorics

import MutableArithmetics as MA

using MultivariatePolynomials
const MP = MultivariatePolynomials

function div_number_test()
    Mod.@polyvar x
    @test div(2x, 2) == x
    @test div(6x + 9x^2, 3) == 2x + 3x^2
    @test div(6x + 9x^2, 4) == x + 2x^2
    if VERSION >= v"1.6"
        @test div(10x^3, 4, RoundUp) == 3x^3
        @test div(6x + 9x^2, 4, RoundUp) == 2x + 3x^2
    end
end

function gcd_lcm_test()
    Mod.@polyvar x y z
    @test gcd(x * y, y) == y
    @test lcm(x * y, y) == x * y
    @test gcd(x * y^2 * z, y * z^3) == y * z
    @test lcm(x * y^2 * z, y * z^3) == x * y^2 * z^3
    @test gcd(x^2 * y^7 * z^3, x^4 * y^5 * z^2) == x^2 * y^5 * z^2
    @test lcm(x^2 * y^7 * z^3, x^4 * y^5 * z^2) == x^4 * y^7 * z^3
end

# Taken from
# Ideals, Varieties, and Algorithms
# Cox, Little and O'Shea, Fourth edition
# They have been adapted to the grlex ordering
function divides_test()
    Mod.@polyvar x y z
    @test divides(leading_monomial(x + y), x) # leading_monomial(x+y) will be x^1*y^0 -> tricky test !
    @test !divides(leading_monomial(x^2 + y), x)
    @test divides(x * y, x^2 * y)
    @test divides(x * y, x * y^2)
    @test divides(y * z, x * y * z)
    @test !divides(y * z, x * z)
    @test !divides(x^2 * y, x * y^2)
end

function leading_term_test()
    Mod.@polyvar x y z
    @test @inferred(monic(2x^2 + 4y + 2)) == x^2 + 2y + 1
    # See page 60
    p = 4x * y^2 * z + 4z^2 + 7x^2 * z^2 - 5x^3
    @test leading_coefficient(p) == 7
    @test leading_monomial(p) == x^2 * z^2
    @test leading_term(p) == 7x^2 * z^2
end

function divrem_test()
    Mod.@polyvar x y
    p = x * y^2 + 1
    q = x * y + 1
    @test typeof(div(p, q)) == MA.promote_operation(div, typeof(p), typeof(q))
    @test typeof(rem(p, q)) == MA.promote_operation(rem, typeof(p), typeof(q))
    @test coefficient_type(div(p, q)) == Rational{Int}
    @test coefficient_type(rem(p, q)) == Rational{Int}
    @test (@inferred div(p, q)) == y
    @test (@inferred rem(p, q)) == -y + 1
    @test (@inferred div(x * y^2 + x, y)) == x * y
    @test (@inferred rem(x * y^2 + x, y)) == x
    @test (@inferred rem(x^4 + x^3 + (1 + 1e-10) * x^2 + 1, x^2 + x + 1)) == 1
    @test (@inferred rem(
        x^4 + x^3 + (1 + 1e-10) * x^2 + 1,
        x^2 + x + 1;
        ztol = 1e-11,
    )) ≈ -((1 + 1e-10) - 1)x + 1
end

function testdiv(p, ps)
    q, r = @inferred divrem(p, ps)
    @test p == dot(q, ps) + r
    return q, r
end

function multi_div_test()
    Mod.@polyvar x y
    # Example 1
    q, r = testdiv(x * y^2 + 1, [x * y + 1, y + 1])
    @test q == [y, -1]
    @test r == 2
    # Example 2
    q, r = testdiv(x^2 * y + x * y^2 + y^2, [x * y - 1, y^2 - 1])
    @test q == [x + y, 1]
    @test r == x + y + 1
    # Example 4
    q, r = testdiv(x^2 * y + x * y^2 + y^2, [y^2 - 1, x * y - 1])
    @test q == [x + 1, x]
    @test r == 2x + 1
    # Example 5
    q, r = testdiv(x * y^2 - x, [x * y - 1, y^2 - 1])
    @test q == [y, 0]
    @test r == -x + y
    q, r = testdiv(x * y^2 - x, [y^2 - 1, x * y - 1])
    @test q == [x, 0]
    @test r == 0

    @testset "Issue DynamicPolynomials#62" begin
        p = x^2 + x + 1
        q = rem(p, [x^2 - y])
        @test q == x + y + 1
    end

    @test (@inferred rem(
        x^2 * y + (1 + 1e-10) * x * y + 1,
        [x^2 + x, y + 1],
    )) == 1
    @test (@inferred rem(
        x^2 * y + (1 + 1e-10) * x * y + 1,
        [x^2 + x, y + 1];
        ztol = 1e-11,
    )) == -((1 + 1e-10) - 1)x + 1
end

function lcm_test()
    Mod.@polyvar x
    l = @inferred lcm(x^2 - 1, x^2 + 2x + 1)
    exp = -(x - 1) * (x + 1)^2
    @test l == exp || l == -exp
    l = @inferred lcm(x^2 + 2x + 1, x^2 - 1)
    @test l == exp || l == -exp
    @test x^2 - 1 == @inferred lcm(x^2 - 1, x - 1)
    @test x^2 - 1 == @inferred lcm(x - 1, x^2 - 1)
    @test 0 == @inferred lcm(x - x, x + 1)
    @test 0 == @inferred lcm(x + 1, x - x)
    @test 0 == @inferred lcm(x - x, x^2 - x^2)
    @test 0 == @inferred lcm(x^2 - x^2, x - x)
end

@testset "Division" begin
    @testset "div by number" begin
        div_number_test()
    end
    @testset "GCD and LCM" begin
        gcd_lcm_test()
    end
    @testset "Divides" begin
        divides_test()
    end
    @testset "Leading function" begin
        leading_term_test()
    end
    @testset "Division examples" begin
        divrem_test()
    end
    @testset "Division by multiple polynomials examples" begin
        multi_div_test()
    end
    @testset "lcm" begin
        lcm_test()
    end
end

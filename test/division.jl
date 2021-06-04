using LinearAlgebra, Test
using Combinatorics

using MultivariatePolynomials
const MP = MultivariatePolynomials

function gcd_lcm_test()
    Mod.@polyvar x y z
    @test gcd(x*y, y) == y
    @test lcm(x*y, y) == x*y
    @test gcd(x*y^2*z, y*z^3) == y*z
    @test lcm(x*y^2*z, y*z^3) == x*y^2*z^3
    @test gcd(x^2*y^7*z^3, x^4*y^5*z^2) == x^2*y^5*z^2
    @test lcm(x^2*y^7*z^3, x^4*y^5*z^2) == x^4*y^7*z^3
end

# Taken from
# Ideals, Varieties, and Algorithms
# Cox, Little and O'Shea, Fourth edition
# They have been adapted to the grlex ordering
function divides_test()
    Mod.@polyvar x y z
    @test divides(leadingmonomial(x+y), x) # leadingmonomial(x+y) will be x^1*y^0 -> tricky test !
    @test !divides(leadingmonomial(x^2+y), x)
    @test divides(x*y, x^2*y)
    @test divides(x*y, x*y^2)
    @test divides(y*z, x*y*z)
    @test !divides(y*z, x*z)
    @test !divides(x^2*y, x*y^2)
end

function leading_term_test()
    Mod.@polyvar x y z
    @test @inferred(monic(2x^2 + 4y + 2)) == x^2 + 2y + 1
    # See page 60
    p = 4x*y^2*z + 4z^2 + 7x^2*z^2 - 5x^3
    @test leadingcoefficient(p) == 7
    @test leadingmonomial(p) == x^2*z^2
    @test leadingterm(p) == 7x^2*z^2
end

function divrem_test()
    Mod.@polyvar x y
    @test (@inferred div(x*y^2 + 1, x*y + 1)) == y
    @test (@inferred rem(x*y^2 + 1, x*y + 1)) == -y + 1
    @test (@inferred div(x*y^2 + x, y)) == x*y
    @test (@inferred rem(x*y^2 + x, y)) == x
    @test (@inferred rem(x^4 + x^3 + (1+1e-10)*x^2 + 1, x^2 + x + 1)) == 1
    @test (@inferred rem(x^4 + x^3 + (1+1e-10)*x^2 + 1, x^2 + x + 1; ztol=1e-11)) ≈ -((1+1e-10)-1)x + 1
end

function testdiv(p, ps)
    q, r = @inferred divrem(p, ps)
    @test p == dot(q, ps) + r
    q, r
end

function multi_div_test()
    Mod.@polyvar x y
    # Example 1
    q, r = testdiv(x*y^2 + 1, [x*y + 1, y + 1])
    @test q == [y, -1]
    @test r == 2
    # Example 2
    q, r = testdiv(x^2*y + x*y^2 + y^2, [x*y - 1, y^2 - 1])
    @test q == [x + y, 1]
    @test r == x + y + 1
    # Example 4
    q, r = testdiv(x^2*y + x*y^2 + y^2, [y^2 - 1, x*y - 1])
    @test q == [x + 1, x]
    @test r == 2x + 1
    # Example 5
    q, r = testdiv(x*y^2 - x, [x*y - 1, y^2 - 1])
    @test q == [y, 0]
    @test r == -x + y
    q, r = testdiv(x*y^2 - x, [y^2 - 1, x*y - 1])
    @test q == [x, 0]
    @test r == 0

    @testset "Issue DynamicPolynomials#62" begin
        p = x^2 + x + 1
        q = rem(p, [x^2-y])
        @test q == x + y + 1
    end

    @test (@inferred rem(x^2*y + (1+1e-10)*x*y + 1, [x^2 + x, y + 1])) == 1
    @test (@inferred rem(x^2*y + (1+1e-10)*x*y + 1, [x^2 + x, y + 1]; ztol=1e-11)) == -((1+1e-10)-1)x + 1
end

function univariate_gcd_test()
    Mod.@polyvar x
    @test -(x + 1) == @inferred gcd(x^2 - 1, x^2 + 2x + 1)
    @test -(x + 1) == @inferred gcd(x^2 + 2x + 1, x^2 - 1)
    @test -(x + 1) == @inferred gcd(x^2 - 1, x + 1)
    @test -(x + 1) == @inferred gcd(x + 1, x^2 - 1)
    @test   x + 1  == @inferred gcd(x - x, x + 1)
    @test   x + 1  == @inferred gcd(x + 1, x - x)
    @test       0  == @inferred gcd(x - x, x^2 - x^2)
    @test       0  == @inferred gcd(x^2 - x^2, x - x)
end

function _mult_test(a::Number, b)
    return @test iszero(maxdegree(b))
end
function _mult_test(a, b)
    @test iszero(rem(a, b))
    @test iszero(rem(b, a))
end
function mult_test(expected, a, b)
    g = @inferred gcd(a, b)
    @test g isa promote_type(polynomialtype(a), polynomialtype(b))
    _mult_test(expected, g)
end
function sym_test(a, b, g)
    mult_test(g, a, b)
    mult_test(g, b, a)
end
function triple_test(a, b, c)
    sym_test(a * c, b * c, gcd(a, b) * c)
    sym_test(b * a, c * a, gcd(b, c) * a)
    sym_test(c * b, a * b, gcd(c, a) * b)
end

function multivariate_gcd_test(::Type{T}) where {T}
    Mod.@polyvar x y z
    o = one(T)
    # Inspired from https://github.com/JuliaAlgebra/MultivariatePolynomials.jl/issues/160
    f1 = o * x * y + o * x
    f2 = o * y^2
    f3 = o * x
    sym_test(f1, f2, 1)
    sym_test(f2, f3, 1)
    sym_test(f3, f1, x)
    triple_test(f1, f2, f3)

    p1 = z - z
    p2 = z
    @test gcd(p1, p2) == z
    p1 = y - y
    p2 = z
    @test gcd(p1, p2) == z
    p1 = 2*y*z - 1
    p2 = y - z
    @test gcd(p1, p2) == -1
    p1 = 2*y^2*z - y
    p2 = y - z
    @test gcd(p1, p2) == -1
    p1 = 2*y^2*z - y
    p2 = y^3*z - y - z
    @test gcd(p1, p2) == -1
    p1 = 2*y^3 + 2*y^2*z - y
    p2 = y^3 + y^3*z - y - z
    g = @inferred gcd(p1, p2)
    @test g == -1 || g == 1
    p1 = z^4 + 2*y^3 + 2*y^2*z - y
    p2 = y*z^4 + y^3 + y^3*z - y - z
    g = @inferred gcd(p1, p2)
    @test g == -1 || g == 1
    p1 = y*z^3 + z^4 + 2*y^3 + 2*y^2*z - y
    p2 = y^2*z^3 + y*z^4 + y^3 + y^3*z - y - z
    g = @inferred gcd(p1, p2)
    @test g == -1 || g == 1
    p1 = -3*y^2*z^3 - 3*y^4 + y*z^3 + z^4 + 2*y^3 + 2*y^2*z - y
    p2 = 3*y^3*z^3 - 2*y^5 + y^2*z^3 + y*z^4 + y^3 + y^3*z - y - z
    g = @inferred gcd(p1, p2)
    @test g == -1 || g == 1
    p1 = -z^6 - 3*y^2*z^3 - 3*y^4 + y*z^3 + z^4 + 2*y^3 + 2*y^2*z - y
    p2 = -y*z^6 - 3*y^3*z^3 - 2*y^5 + y^2*z^3 + y*z^4 + y^3 + y^3*z - y - z
    g = @inferred gcd(p1, p2)
    @test g == -1 || g == 1
    p1 = -z^6 - 3*y^2*z^3 - 3*y^4 + y*z^3 + z^4 + 2*y^3 + 2*y^2*z - y
    p2 = -y^2*z^6 - 3*y^4*z^3 - 2*y^6 + y^3*z^3 + y^2*z^4 + y^5 + y^4*z - y^2 - y*z
    g = @inferred gcd(p1, p2)
    @test g == -1 || g == 1

    a = (o * x + o * y^2) * (o * z^3 + o * y^2 + o * x)
    b = (o * x + o * y + o * z) * (o * x^2 + o * y)
    c = (o * x + o * y + o * z) * (o * z^3 + o * y^2 + o * x)
    sym_test(a, b, 1)
    sym_test(b, c, x + y + z)
    sym_test(c, a, z^3 + y^2 + x)
    triple_test(a, b, c)
end

function lcm_test()
    Mod.@polyvar x
    @test -(x-1) * (x+1)^2 == @inferred lcm(x^2 - 1, x^2 + 2x + 1)
    @test -(x-1) * (x+1)^2 == @inferred lcm(x^2 + 2x + 1, x^2 - 1)
    @test x^2 - 1 == @inferred lcm(x^2 - 1, x - 1)
    @test x^2 - 1 == @inferred lcm(x - 1, x^2 - 1)
    @test       0 == @inferred lcm(x - x, x + 1)
    @test       0 == @inferred lcm(x + 1, x - x)
    @test       0 == @inferred lcm(x - x, x^2 - x^2)
    @test       0 == @inferred lcm(x^2 - x^2, x - x)
end

function deflation_test()
    Mod.@polyvar a b c
    function _test(p, eshift, edefl)
        shift, defl = MP.deflation(p)
        @test shift == eshift
        @test defl == edefl
        q = MP.deflate(p, shift, defl)
        @test p == MP.inflate(q, shift, defl)
    end
    for (x, y, z) in permutations((a, b, c))
        _test(x + x^2 + y * x, x, x * y)
        p = x^2 + x^2 * y^3 + y^6 * x^4
        _test(p, x^2, x^2 * y^3)
        @test p == MP.deflate(p, z^0, z^0)
        _test(p, x^2, x^2 * y^3)
    end
end

function extracted_variable_test()
    Mod.@polyvar x y z
    function _test(p1, p2, w1, w2)
        i1, i2, n = MP._extracted_variable(p1, p2)
        v(p, i) = iszero(i) ? nothing : variables(p)[i]
        if string(Mod) == "DynamicPolynomials"
            alloc_test(() -> MP._extracted_variable(p1, p2), 0)
        end
        @test v(p1, i1) == w1
        @test v(p2, i2) == w2
        i2, i1, n = MP._extracted_variable(p2, p1)
        if string(Mod) == "DynamicPolynomials"
            alloc_test(() -> MP._extracted_variable(p2, p1), 0)
        end
        @test v(p1, i1) == w1
        @test v(p2, i2) == w2
    end
    _test(x - x, y - y, nothing, nothing)
    _test(x + y + z, x - z, y, nothing)
    _test(x + y + z, x - y, z, nothing)
    _test(x + y + z, y - z, x, nothing)
    _test(x^4 + y^5 + z^6, x^3 + y^2 + z, z, z)
    _test(x^6 + y^5 + z^4, x + y^2 + z^3, x, x)
    _test(x^6 + y + z^4 - y, x + y + z^3 - y, x, x)
end

@testset "Division" begin
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
    @testset "Univariate gcd" begin
        univariate_gcd_test()
    end
    @testset "Multivariate gcd $T" for T in [Int, Rational{BigInt}]
        multivariate_gcd_test(T)
    end
    @testset "lcm" begin
        lcm_test()
    end
end

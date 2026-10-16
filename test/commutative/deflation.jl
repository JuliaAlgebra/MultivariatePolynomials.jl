using Test
using Combinatorics
using MultivariatePolynomials
import MultivariatePolynomials as MP
if !@isdefined(Mod)
    import DynamicPolynomials as Mod
end

@testset "Deflation" begin
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

    @testset "Different variable sets" begin
        for (x, y) in permutations((a, b))
            m, shift, defl = y^2, y^2, x * y
            @test isone(@inferred MP.deflate(m, shift, defl))
            @test m == shift == y^2
            @test collect(variables(m)) == collect(variables(shift)) == [y]
            @test defl == x * y

            m, shift, defl = y^2, monomial(x), y^2
            @test (@inferred MP.inflate(m, shift, defl)) == x * y^4
            @test m == defl == y^2
            @test collect(variables(m)) == collect(variables(defl)) == [y]
            @test shift == x
            @test collect(variables(shift)) == [x]
        end
    end
end

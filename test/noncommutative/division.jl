using Test
using MultivariatePolynomials
import MutableArithmetics as MA

if !@isdefined(Mod)
    import DynamicPolynomials as Mod
end

@testset "Noncommutative remainder" begin
    Mod.@ncpolyvar x y z
    relation = y^2 - (2x - 1)
    @test rem(2x + 3y^2, relation) == 8x - 3
    @test rem(x * y^3 + y^2, [relation]) == 2x^2 * y - x * y + 2x - 1
    @test rem(x * y^2 * z, relation) == 2x^2 * z - x * z
    @test rem(y^2 * z, relation) == 2x * z - z
    @test rem(x * y^2, relation) == 2x^2 - x
    @test rem(y * x * y, relation) == y * x * y
    @test rem(y^4, relation) == 4x^2 - 4x + 1
    @test rem(y^3, relation) == 2x * y - y
    @test rem(y^5, relation) == 4x^2 * y - 4x * y + y
    @test rem(y^2, 2y^2 - x - 1) == (x + 1) / 2
    @test rem(x * y * z, [x * y - x, y * z - z]) == x * z
    @test rem(x * y * z, [y * z - z, x * y - x]) == x * z
    @test rem(x * y * z, [x * y - x, y * z - y]) == x * z
    @test rem(x * y * z, [y * z - y, x * y - x]) == x
    @test rem(y^2, [relation, x - 1]) == 1
    @test rem(x * y, polynomial(2x^0)) == 0
    @test rem(zero(x), relation) == 0
    @test rem(x + y, typeof(relation)[]) == x + y
    @test rem(1e-12 * x + y, relation; ztol = 1e-10) == y
    @test rem(1e-12 * x + y, relation; ztol = 1e-14) == 1e-12 * x + y
    @test_throws DivideError rem(x, zero(x))
    @test_throws DivideError rem(zero(x), [zero(x)])

    # Padding and repeated adjacent variables represent the same word.
    padded = monomial([x, y, x, z], [0, 1, 0, 1])
    split = monomial([y, y], [1, 1])
    @test rem(padded, y * z - x) == x
    @test rem(split, relation) == 2x - 1
    @test rem(x * y * x * y^2 * z, relation) == 2x * y * x^2 * z - x * y * x * z

    for T in (Int, BigInt, Rational{Int}, Float64)
        f = T(3) * x * y^2 * z + T(2) * z
        g = [T(2) * y^2 - T(1) * x - T(1)]
        f_copy, g_copy = MA.mutable_copy(f), MA.mutable_copy(g)
        result = rem(f, g)
        @test result == (3 // 2) * x^2 * z + (3 // 2) * x * z + 2z
        @test rem(result, g) == result
        @test f == f_copy
        @test g == g_copy
    end
end

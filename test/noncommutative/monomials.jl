using Test
using MultivariatePolynomials
if !@isdefined(Mod)
    import DynamicPolynomials as Mod
end

function monomial_word(m)
    return Tuple(v for (v, e) in powers(m) for _ in 1:e)
end

@testset "Noncommutative monomial enumeration" begin
    Mod.@ncpolyvar x y z
    @test (@inferred monomials((x, y), 2)) == [y * x, y^2, x * y, x^2]
    @test (@inferred monomials([x, y], 0:2)) ==
          [x^0, y, x, y * x, y^2, x * y, x^2]
    cubes = @inferred monomials([x, y], 3)
    @test cubes ==
          [y * x * y, y * x^2, y^2 * x, y^3, x * y * x, x * y^2, x^2 * y, x^3]
    @test collect.(exponents.(cubes)) == reverse([
        [3, 0, 0, 0],
        [2, 1, 0, 0],
        [1, 2, 0, 0],
        [1, 1, 1, 0],
        [0, 3, 0, 0],
        [0, 2, 1, 0],
        [0, 1, 2, 0],
        [0, 1, 1, 1],
    ])
    @test all(m -> variables(m) == [x, y, x, y], cubes)

    # Compare the generated words to the Cartesian product, independently of
    # the padded exponent representation used by the generator.
    for vars in ([x], [x, y], [x, y, z]), deg in 0:4
        monos = @inferred monomials(vars, deg)
        expected = Set(Iterators.product(ntuple(_ -> vars, deg)...))
        @test length(monos) == length(vars)^deg
        @test Set(monomial_word.(monos)) == expected
        @test issorted(monos)
    end

    vars, degrees = [y, x, y, x, y], [3, 0, 1, 3]
    monos = @inferred monomials(vars, degrees)
    @test length(monos) == 11
    @test issorted(monos)
    @test monos == monomials([x, y], [0, 1, 3])
    @test (@inferred monomials(@view(vars[1:2]), degrees)) == monos
    @test vars == [y, x, y, x, y]
    @test degrees == [3, 0, 1, 3]
    @test (@inferred monomials((x, y), 2, m -> m == y * x)) == [y * x]
    @test isempty(@inferred monomials(vars, 0:2, _ -> false))
    @test isempty(@inferred monomials(vars, Int[]))
    @test isempty(@inferred monomials(vars, 1:0))
    @test (@inferred monomials(x, 0:3)) == [x^0, x, x^2, x^3]
    @test (@inferred monomials(typeof(x)[], 0:2)) == [x^0]
    @test isempty(@inferred monomials(typeof(x)[], 1:2))
    @test_throws ArgumentError monomials(vars, -1)
    @test_throws ArgumentError monomials(vars, -1:1)
    @test_throws ArgumentError monomials(vars, [1, -1])

    seen = typeof(x^0)[]
    keep_yx = m -> begin
        push!(seen, m)
        return m == y * x
    end
    kept = @inferred monomials([x, y], 2, keep_yx)
    @test seen == [y * x, y^2, x * y, x^2]
    @test kept == [y * x]
    exponents(seen[1])[1] = 5
    @test seen[2:4] == [y^2, x * y, x^2]
    @test variables(x) == [x]
    @test variables(y) == [y]
end

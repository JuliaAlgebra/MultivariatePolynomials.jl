import Test: @inferred

@testset "Substitution" begin
    Mod.@polyvar x[1:3]

    @test subs(2, x[1]=>3) == 2
    @test iszero((x[1]-x[1])(x[1]=>x[2]))
    @test subs(CustomPoly(x[1]+x[2]), x[2]=>x[1]) == 2x[1]

    a = (x[1])(x[1]=>x[2])
    b = x[2]
    @test (x[1])(x[1]=>x[2]) == x[2]

    p = x[1]*x[2]*x[3]
    @test convert(Int, p(x => (1, 2, 3))) == 6

    p = x[1]^2 + x[1]*x[3] - 3
    @test convert(Int, p(x => (5, x[1]*x[2], 4))) == 42

    p = x[1]^2 + x[2]^2
    q = p(x[1:2] => [1 -1; 1 1] * vec(x[1:2]))
    @test q == 2p

    t = 2.0 * x[1] * x[2]^2
    @test t(x => [1, 2, 3]) == 8.0
    @test t((x[2], x[1]) => [1, 3]) == 6.0

    p = x[1] + x[2] + 2*x[1]^2 + 3*x[1]*x[2]^2
    #@inferred p((x[1], x[2]) => (1.0, 2.0))

    @inferred subs(p, x[2] => 2.0)
    @test subs(p, x[2] => 2.0) == 13x[1] + 2 + 2x[1]^2
    @inferred subs(p, x[2] => x[1])
    @inferred subs(p, x[2] => x[1]) == 2x[1] + 2x[1]^2 + 3x[1]^3

    p = x[1] - x[1]
    @test iszero(p)
    @inferred p(x[1] => 1)
    @test iszero(p(x[1] => 1))
    @inferred subs(p, x[1] => 1)
    @test iszero(subs(p, x[1] => 1))

    q = (x[1] + 1) / (x[1] + 2)
    @test isapproxzero(q(x[1] => -1))
    @test !isapproxzero(q(x[1] => 1))
    @test q(x[1] => 1) ≈ 2/3

    q = (x[1] + x[3]) / (x[2] - 1)
    @test subs(q, x[1:2] => (x[2], x[1])) == (x[2] + x[3]) / (x[1] - 1)

    # Taken from perf/runbenchmark.jl
    p = x[1] + x[2] + 2x[1]^2 + x[2]^3

    varst = (x[1], x[2])
    varsv = [x[1], x[2]]
    valst = (1.0, 2.0)
    valsv = [1.0, 2.0]

    @test variable.(Tuple(varsv)) == varst
    @test variable.(Tuple(varsv)) isa Tuple{<:AbstractVariable, <:AbstractVariable}

    @test p(variables(p) => valst) == 13
    @test p(variables(p) => valsv) == 13
    @test p(varst => valst) == 13
    @test p(varsv => valsv) == 13
end

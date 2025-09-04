module TestComparison

using Test
using MultivariatePolynomials

function _test(object, M; kws...)
    it = ExponentsIterator{M}(object; kws...)
    v = collect(Iterators.take(it, 20))
    @test issorted(v, lt = (a, b) -> cmp(M(), a, b) < 0)
end

function _test(nvars::Int, M; kws...)
    _test(zeros(Int, nvars), M; inline = false, kws...)
    _test(zeros(Int, nvars), M; inline = true, kws...)
    _test(ntuple(zero, nvars), M; inline = false, kws...)
    return
end

function test_errors()
    err = ArgumentError(
        "The `mindegree` of `ExponentsIterator` cannot be negative.",
    )
    @test_throws err ExponentsIterator{LexOrder}([0], mindegree = -1)
    M = Reverse{LexOrder}
    err = ArgumentError(
        "Ordering `$M` is not a valid ordering, use `Graded{$M}` instead.",
    )
    @test_throws err ExponentsIterator{M}([0], maxdegree = 2)
    exps = ExponentsIterator{LexOrder}([0])
    err = ErrorException(
        "The iterator is infinity because `maxdegree` is `nothing`.",
    )
    @test_throws err length(exps)
end

function test_exponents_iterator()
    @testset "nvariables = $nvars" for nvars in 0:3
        @testset "mindegree = $mindegree" for mindegree in 0:3
            @testset "maxdegree = $maxdegree" for maxdegree in
                                                  vcat(nothing, -1:3)
                for L in [LexOrder, InverseLexOrder]
                    @testset "M = $M" for M in
                                          [L, Graded{L}, Graded{Reverse{L}}]
                        _test(nvars, M; mindegree, maxdegree)
                    end
                end
            end
        end
    end
end

function test_compare()
    lex = LexOrder()
    grlex = Graded{LexOrder}()
    rinvlex = Reverse{InverseLexOrder}()
    grevlex = Graded{Reverse{InverseLexOrder}}()
    @test cmp(grlex, [1, 0, 1], [1, 1, 0]) == -1
    @test cmp(grlex, [1, 1, 0], [1, 0, 1]) == 1
    # [CLO13, p. 58]
    @test cmp(lex, 1:3, [3, 2, 0]) < 0
    @test cmp(grlex, 1:3, [3, 2, 0]) > 0
    @test cmp(rinvlex, 1:3, [3, 2, 0]) < 0
    @test cmp(grevlex, 1:3, [3, 2, 0]) > 0
    @test cmp(lex, [1, 2, 4], [1, 1, 5]) > 0
    @test cmp(grlex, [1, 2, 4], [1, 1, 5]) > 0
    @test cmp(rinvlex, [1, 2, 4], [1, 1, 5]) > 0
    @test cmp(grevlex, [1, 2, 4], [1, 1, 5]) > 0
    # [CLO13, p. 59]
    @test cmp(lex, (5, 1, 1), (4, 1, 2)) > 0
    @test cmp(grlex, (5, 1, 1), (4, 1, 2)) > 0
    @test cmp(rinvlex, (5, 1, 1), (4, 1, 2)) > 0
    @test cmp(grevlex, (5, 1, 1), (4, 1, 2)) > 0
    # [CLO13] Cox, D., Little, J., & OShea, D.
    # *Ideals, varieties, and algorithms: an introduction to computational algebraic geometry and commutative algebra*.
    # Springer Science & Business Media, **2013**.
end

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

end # module

TestComparison.runtests()

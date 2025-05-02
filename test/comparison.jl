module TestComparison

using Test
using MultivariatePolynomials

function _test(object, M; kws...)
    it = ExponentsIterator{M}(object; kws...)
    v = collect(Iterators.take(it, 20))
    @test issorted(v, lt = (a, b) -> compare(a, b, M) < 0)
end

function _test(nvars::Int, M; kws...)
    _test(zeros(Int, nvars), M; inline = false, kws...)
    _test(zeros(Int, nvars), M; inline = true, kws...)
    _test(ntuple(zero, nvars), M; inline = false, kws...)
    return
end

function test_exponents_iterator()
    @testset "nvariables = $nvars" for nvars in 0:3
        @testset "mindegree = $mindegree" for mindegree in 0:3
            @testset "maxdegree = $maxdegree" for maxdegree in
                                                  vcat(nothing, 0:3)
                for L in [LexOrder, InverseLexOrder]
                    @testset "M = $M" for M in [L, Graded{L}]
                        _test(nvars, M; mindegree, maxdegree)
                    end
                end
            end
        end
    end
end

function test_compare()
    lex = LexOrder
    grlex = Graded{lex}
    rinvlex = Reverse{InverseLexOrder}
    grevlex = Graded{rinvlex}
    @test compare([1, 0, 1], [1, 1, 0], grlex) == -1
    @test compare([1, 1, 0], [1, 0, 1], grlex) == 1
    # [CLO13, p. 58]
    @test compare(1:3, [3, 2, 0], lex) < 0
    @test compare(1:3, [3, 2, 0], grlex) > 0
    @test compare(1:3, [3, 2, 0], rinvlex) < 0
    @test compare(1:3, [3, 2, 0], grevlex) > 0
    @test compare([1, 2, 4], [1, 1, 5], lex) > 0
    @test compare([1, 2, 4], [1, 1, 5], grlex) > 0
    @test compare([1, 2, 4], [1, 1, 5], rinvlex) > 0
    @test compare([1, 2, 4], [1, 1, 5], grevlex) > 0
    # [CLO13, p. 59]
    @test compare((5, 1, 1), (4, 1, 2), lex) > 0
    @test compare((5, 1, 1), (4, 1, 2), grlex) > 0
    @test compare((5, 1, 1), (4, 1, 2), rinvlex) > 0
    @test compare((5, 1, 1), (4, 1, 2), grevlex) > 0
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

module TestExponentsIterator

using Test
using MultivariatePolynomials

function _test(object, M; kws...)
    it = ExponentsIterator{M}(object; kws...)
    v = collect(Iterators.take(it, 20))
    @test issorted(v, lt = (a, b) -> compare(a, b, M))
end

function _test(nvars::Int, M; kws...)
    _test(zeros(Int, nvars); inline = false, kws...)
    _test(zeros(Int, nvars); inline = true, kws...)
    _test(ntuple(zero, nvars); inline = false, kws...)
    return
end

function test_all()
    @testset "nvariables = $nvars" for nvars in 0:3
        @testset "mindegree = $mindegree" for mindegree in 0:3
            @testset "maxdegree = $maxdegree" for maxdegree in vcat(nothing; 0:3)
                for L in [LexOrder, InverseLexOrder]
                    @testset "M = $M" for M in [
                        L, Graded{L},
                    ]
                        _test(nvars, M; mindegree, maxdegree)
                    end
                end
            end
        end
    end
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

TestExponentsIterator.runtests()

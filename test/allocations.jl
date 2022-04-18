module TestAllocations

include("utils.jl")

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
end

using MutableArithmetics
using TypedPolynomials

function test_polynomial_merge()
    @polyvar x
    p = x^2 + x + 1
    for q in [x, x^2, 2x, 2x + 3]
        alloc_test(0) do
            add!!(p, q)
        end
        alloc_test(0) do
            sub!!(p, q)
        end
    end
end

end
TestAllocations.runtests()

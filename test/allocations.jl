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
import MultivariatePolynomials
const MP = MultivariatePolynomials
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


function test_promotion()
    @polyvar x y
    p = x * y + x + y + 1
    alloc_test(0) do
        promote_operation(MP.substitute, MP.Subs, typeof(p), Pair{typeof(x),Int})
    end

function test_isapproxzero()
    @polyvar x y
    p = x * y + x + y + 1
    alloc_test(0) do
        isapproxzero(p)
    end
    q = 1e-10 * x * y + 1e-12 * x
    alloc_test(0) do
        isapproxzero(q)
    end
    alloc_test(0) do
        isapproxzero(q; ztol = 1e-8)
    end
end

function _test_gcd(T)
    o = one(T)
    @polyvar x y
    t1 = o * x * y
    t2 = o * x
    alloc_test(0) do
        gcd(t1, t2)
    end
end

function test_gcd()
    _test_gcd(Int)

end

end
TestAllocations.runtests()

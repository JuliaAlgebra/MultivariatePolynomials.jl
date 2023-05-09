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
    # Want to avoid it converting everything to polynomial
    Q = AbstractPolynomialLike[x, x^2, 2x, 2x+3]
    for q in Q
        alloc_test(0) do
            add!!(p, q)
        end
        alloc_test(0) do
            sub!!(p, q)
        end
        for r in Q
            if q isa AbstractPolynomial && r isa AbstractPolynomial
                continue
            end
            alloc_test(0) do
                add_mul!!(p, q, r)
            end
            alloc_test(0) do
                sub_mul!!(p, q, r)
            end
        end
    end
end

function __test_leading(p)
    alloc_test(0) do
        leading_term(p)
    end
    alloc_test(0) do
        leading_coefficient(p)
    end
    alloc_test(0) do
        leading_monomial(p)
    end
end
function _test_leading(T)
    @polyvar x y
    __test_leading(T(4) * x^2 * y + x * y + 2x)
    __test_leading(T(4) * x^2 * y)
end

function test_leading()
    _test_leading(Int)
    _test_leading(BigInt)
end

function test_promotion()
    @polyvar x y
    p = x * y + x + y + 1
    alloc_test(0) do
        promote_operation(MP.substitute, MP.Subs, typeof(p), Pair{typeof(x),Int})
    end
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
    alloc_test(0) do
        iszero(q)
    end
    alloc_test(0) do
        isone(q)
    end
end

function _test_term_gcd(T)
    o = one(T)
    @polyvar x y
    t1 = o * x * y
    t2 = o * x
    alloc_test(0) do
        gcd(t1, t2)
    end
end

function test_term_gcd()
    _test_term_gcd(Int)
end

function _pseudo_rem_test(p1, p2, algo)
    backup_1 = deepcopy(p1)
    MA_copy_1 = mutable_copy(p1)
    backup_2 = deepcopy(p2)
    MA_copy_2 = mutable_copy(p2)
    buffer = buffer_for(MP.rem_or_pseudo_rem, typeof(p1), typeof(p2), typeof(algo))
    mutable_alloc_test(p1, 0) do p1
        buffered_operate!!(buffer, MP.rem_or_pseudo_rem, p1, p2, algo)
    end
    @test backup_1 == MA_copy_1
    @test backup_2 == MA_copy_2
    @test backup_2 == p2
    return MA_copy_1, p2
end

function _test_div(T)
    o = one(T)
    @polyvar x y
    p1 = 1o * x^2 + 3o * x + 1o
    p2 = 1o * x + 1o
    @testset "primitive_rem=$primitive_rem" for primitive_rem in [false, true]
        @testset "skip_last=$skip_last" for skip_last in [false, true]
            algo = GeneralizedEuclideanAlgorithm(primitive_rem, skip_last)
            p1, p2 = _pseudo_rem_test(p1, p2, algo)
        end
    end
end

function test_div()
    @testset "$T" for T in [Int] # FIXME , BigInt] # TODO add Rational{BigInt} when https://github.com/jump-dev/MutableArithmetics.jl/issues/187 is fixed
        _test_div(T)
    end
end

function _test_univ_gcd(T, algo)
    o = one(T)
    @polyvar x
    p1 = o * x^2 + 2o * x + o
    p2 = o * x + o
    mutable_alloc_test(mutable_copy(p1), 0) do p1
        MP.univariate_gcd(p1, p2, algo, IsMutable(), IsMutable())
    end
end

function test_univ_gcd()
    algo = GeneralizedEuclideanAlgorithm()
    _test_univ_gcd(Int, algo)
    #_test_univ_gcd(BigInt, algo) # TODO
end

end

TestAllocations.runtests()

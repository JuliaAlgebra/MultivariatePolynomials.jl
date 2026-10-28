if !@isdefined(Mod)
    import DynamicPolynomials as Mod
    include("../utils.jl")
end

using LinearAlgebra, Test
using MultivariatePolynomials
using ChainRulesCore

function test_chain_rule(dot, op, args, Δin, Δout)
    output = op(args...)
    foutput, fΔout = ChainRulesCore.frule((NoTangent(), Δin...), op, args...)
    @test output == foutput
    routput, pullback = ChainRulesCore.rrule(op, args...)
    @test output == routput
    rΔin = pullback(Δout)
    @test rΔin[1] == NoTangent()
    @test dot(Δin, rΔin[2:end]) ≈ dot(fΔout, Δout)
end

function coefficient_dot(p::AbstractPolynomialLike, q::AbstractPolynomialLike)
    return sum(terms(p); init = 0.0im) do t
        return conj(coefficient(t)) * coefficient(q, monomial(t))
    end
end
coefficient_dot(p::Tuple, q::Tuple) = sum(map(coefficient_dot, p, q))

@testset "ChainRulesCore" begin
    Mod.@polyvar x y
    p = 1.1x + y
    q = (-0.1 + im) * x - y

    output, pullback = ChainRulesCore.rrule(+, q)
    @test output == q
    @test pullback(2) == (NoTangent(), 2)
    @test pullback(x + 3) == (NoTangent(), x + 3)

    output, pullback = ChainRulesCore.rrule(-, q)
    @test output ≈ -q
    @test pullback(2) == (NoTangent(), -2)
    @test pullback(x + 3im) == (NoTangent(), -(x + 3im))

    output, pullback = ChainRulesCore.rrule(+, p, q)
    @test output == (1.0 + im)x
    @test pullback(2) == (NoTangent(), 2, 2)
    @test pullback(x + 3) == (NoTangent(), x + 3, x + 3)

    output, pullback = ChainRulesCore.rrule(-, p, q)
    @test output ≈ (1.2 - im) * x + 2y
    @test pullback(2) == (NoTangent(), 2, -2)
    @test pullback(1im * x + 3) == (NoTangent(), 1im * x + 3, -1im * x - 3)

    # Multiplication is a convolution of coefficients; its adjoint is a
    # correlation. Multiplying the output cotangent by q' would give x^3.
    output, pullback =
        ChainRulesCore.rrule(*, polynomial(1.0x), polynomial(1.0x))
    @test output == x^2
    @test pullback(1.0x^2) == (NoTangent(), x, x)
    @test pullback(1.0y) == (NoTangent(), 0, 0)
    @test pullback(ZeroTangent()) == (NoTangent(), 0, 0)
    @test pullback(Thunk(() -> 1.0x^2)) == (NoTangent(), x, x)

    output, pullback = ChainRulesCore.rrule(differentiate, p, x)
    @test output == 1.1
    @test pullback(q) == (NoTangent(), (-0.2 + 2im) * x^2 - x * y, NoTangent())
    @test pullback(1x) == (NoTangent(), 2x^2, NoTangent())

    test_chain_rule(dot, +, (p,), (q,), p)
    test_chain_rule(dot, +, (q,), (p,), q)

    test_chain_rule(dot, -, (p,), (q,), p)
    test_chain_rule(dot, -, (p,), (p,), q)

    test_chain_rule(dot, +, (p, q), (q, p), p)
    test_chain_rule(dot, +, (p, q), (p, q), q)

    test_chain_rule(dot, -, (p, q), (q, p), p)
    test_chain_rule(dot, -, (p, q), (p, q), q)

    test_chain_rule(coefficient_dot, *, (p, q), (q, p), p * q)
    test_chain_rule(coefficient_dot, *, (p, q), (p, q), q * q)
    test_chain_rule(coefficient_dot, *, (q, p), (p, q), q * q)
    test_chain_rule(coefficient_dot, *, (p, q), (q, p), q * q)

    _dot(p, q) = coefficient_dot(p, q)
    function _dot(
        px::Tuple{<:AbstractPolynomial,NoTangent},
        qx::Tuple{<:AbstractPolynomial,NoTangent},
    )
        return _dot(px[1], qx[1])
    end
    test_chain_rule(_dot, differentiate, (p, x), (q, NoTangent()), p)
    test_chain_rule(
        _dot,
        differentiate,
        (p, x),
        (q, NoTangent()),
        differentiate(p, x),
    )
    test_chain_rule(
        _dot,
        differentiate,
        (p, x),
        (q, NoTangent()),
        differentiate(q, x),
    )
    test_chain_rule(_dot, differentiate, (p, x), (p * q, NoTangent()), p)
end

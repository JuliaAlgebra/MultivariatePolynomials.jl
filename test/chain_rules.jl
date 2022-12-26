using LinearAlgebra, Test
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

function _dot(p, q)
    monos = monovec([monomials(p); monomials(q)])
    return dot(coefficient.(p, monos), coefficient.(q, monos))
end
function _dot(px::Tuple, qx::Tuple)
    return _dot(first(px), first(qx)) + _dot(Base.tail(px), Base.tail(qx))
end
function _dot(::Tuple{}, ::Tuple{})
    return MultivariatePolynomials.MA.Zero()
end
function _dot(::NoTangent, ::NoTangent)
    return MultivariatePolynomials.MA.Zero()
end

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
    @test pullback(x + 3im) == (NoTangent(), -x - 3im)

    output, pullback = ChainRulesCore.rrule(+, p, q)
    @test output == (1.0 + im)x
    @test pullback(2) == (NoTangent(), 2, 2)
    @test pullback(x + 3) == (NoTangent(), x + 3, x + 3)

    output, pullback = ChainRulesCore.rrule(-, p, q)
    @test output ≈ (1.2 - im) * x + 2y
    @test pullback(2) == (NoTangent(), 2, -2)
    @test pullback(im * x + 3) == (NoTangent(), im * x + 3, -im * x - 3)

    output, pullback = ChainRulesCore.rrule(differentiate, p, x)
    @test output == 1.1
    @test pullback(q) == (NoTangent(), (-0.2 + 2im) * x^2 - x*y, NoTangent())
    @test pullback(1x) == (NoTangent(), 2x^2, NoTangent())

    for d in [dot, _dot]
        test_chain_rule(d, +, (p,), (q,), p)
        test_chain_rule(d, +, (q,), (p,), q)

        test_chain_rule(d, -, (p,), (q,), p)
        test_chain_rule(d, -, (p,), (p,), q)

        test_chain_rule(d, +, (p, q), (q, p), p)
        test_chain_rule(d, +, (p, q), (p, q), q)

        test_chain_rule(d, -, (p, q), (q, p), p)
        test_chain_rule(d, -, (p, q), (p, q), q)
    end

    test_chain_rule(_dot, *, (p, q), (q, p), p * q)
    test_chain_rule(_dot, *, (p, q), (p, q), q * q)
    test_chain_rule(_dot, *, (q, p), (p, q), q * q)
    test_chain_rule(_dot, *, (p, q), (q, p), q * q)

    test_chain_rule(_dot, differentiate, (p, x), (q, NoTangent()), p)
    test_chain_rule(_dot, differentiate, (p, x), (q, NoTangent()), differentiate(p, x))
    test_chain_rule(_dot, differentiate, (p, x), (q, NoTangent()), differentiate(q, x))
    test_chain_rule(_dot, differentiate, (p, x), (p * q, NoTangent()), p)
end

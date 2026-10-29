if !@isdefined(Mod)
    include("../SimplePolynomials.jl")
    import .SimplePolynomials as Mod
    include("../utils.jl")
end

using LinearAlgebra, Test
using MultivariatePolynomials
using ChainRulesCore
import MutableArithmetics as MA

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
    T = MA.promote_operation(
        *,
        MA.promote_operation(conj, coefficient_type(p)),
        coefficient_type(q),
    )
    return sum(terms(p); init = zero(T)) do t
        return conj(coefficient(t)) * coefficient(q, monomial(t))
    end
end
coefficient_dot(p::Tuple, q::Tuple) = sum(map(coefficient_dot, p, q))

@testset "ChainRulesCore" begin
    Mod.@polyvar x y
    p = 1.1x + y
    q = (-0.1 + im) * x - y

    @testset "Coefficient pairing" begin
        ext = Base.get_extension(
            MultivariatePolynomials,
            :MultivariatePolynomialsChainRulesCoreExt,
        )
        p32, q32 = polynomial(2.0f0 * x), polynomial(3.0f0 * y)
        @test dot(p32, q32) == 6.0f0 * x * y
        @test ext.coefficient_dot(p32, q32) === 0.0f0
        @test dot(p32, p32) == 4.0f0 * x^2
        @test ext.coefficient_dot(p32, p32) === 4.0f0
        for inner in (coefficient_dot, ext.coefficient_dot)
            p32 = polynomial(2.0f0 * x)
            q32 = polynomial(3.0f0 * y)
            # The stored coefficient vectors have equal lengths, but their
            # entries refer to different monomials.
            @test inner(p32, q32) === 0.0f0
            @test inner(zero(p32), q32) === 0.0f0
            @test inner(p32, zero(q32)) === 0.0f0
            @test inner(p32, p32) === 4.0f0
            @test inner(p32 + q32, q32) === 9.0f0
            @test inner(q32, p32 + q32) === 9.0f0
            @test inner((1.0f0 + 2.0f0im) * x, (3.0f0 + 4.0f0im) * x) ===
                  11.0f0 - 2.0f0im
        end
    end

    @testset "Symbolic powers" begin
        for m in (x, x * y), n in (0, 1, 3)
            output, pullback = ChainRulesCore.rrule(^, m, n)
            @test output == m^n
            for delta in (1.0, x + y, ZeroTangent(), Thunk(() -> x))
                @test pullback(delta) == (NoTangent(), NoTangent(), NoTangent())
            end
        end
        # Terms still have a differentiable coefficient.
        _, pullback = ChainRulesCore.rrule(^, 2.0x, 2)
        @test pullback(x^2) == (NoTangent(), 4x, NoTangent())
    end

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

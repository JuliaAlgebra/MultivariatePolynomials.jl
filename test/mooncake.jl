module TestMooncake

using Test, Random
using JET, AllocCheck
if !isdefined(Main, :Mod)
    include("SimplePolynomials.jl")
end
const Mod = isdefined(Main, :Mod) ? Main.Mod : SimplePolynomials
import MultivariatePolynomials as MP
import Mooncake as MC
import ChainRulesCore as CRC

const Ext = Base.get_extension(MP, :MultivariatePolynomialsMooncakeExt)

function gradient(f, args...)
    cache = MC.prepare_gradient_cache(f, args...)
    return MC.value_and_gradient!!(cache, f, args...)
end

term_product(a, b, x) = MP.coefficient((a * x) * (b * x), x^2)
term_sum(a, b, x) = MP.coefficient(a * x + b * x, x)
poly_coefficient(c, x) = MP.coefficient(MP.polynomial(c, [one(x), x, x^2]), x)
view_coefficient(c, x) = poly_coefficient(view(c, 1:3), x)
poly_derivative(a, x) = MP.coefficient(MP.differentiate(a * x^2, x), x)
poly_product(c, x) = MP.coefficient((c[1] * x + c[2]) * (c[3] * x + c[4]), x)
poly_power(a, x, n) = MP.coefficient((a * x + 1)^n, x)
term_power(a, x, n) = MP.coefficient((a * x)^n, x^n)
poly_scalar(a, b, x) = MP.coefficient(((a * x + b) / b + b) * b, x)
left_scalar(a, x) = MP.coefficient(a - (2.0 * x + a), x)
integer_constant(a, x) = MP.coefficient(2x, x) + a
term_constructor(a, x) = MP.coefficient(MP.term(a, x))
poly_constructor(a, x) = MP.coefficient(MP.polynomial(MP.term(a, x)), x)
unary(a, x) = MP.coefficient(-(+MP.polynomial(a * x)), x)
function complex_product(c, x)
    return real(MP.coefficient((c[1] * x + c[2]) * (c[3] * x + c[4]), x))
end
function complex_mixed(a, x)
    return real(MP.coefficient((a * x + 1) * ((1.0 + 2im) * x + 3im), x))
end
scaled_coefficient(c, p, x) = real(MP.coefficient(c * p, x))
missing_coefficient(p, x) = MP.coefficient(p, x^3)
product_different_parents(p, q, x, y) = MP.coefficient(p * q, x * y)
product_same_argument(p, q, x) = MP.coefficient(p * q, x^2)
repeated_argument(p, x) = MP.coefficient(p * p + p, x^2)
identity_polynomial(p, x) = MP.coefficient(+MP.polynomial(p), x)

function test_backend(x, y)
    @testset "Scalar and term inputs" begin
        for (a, b) in ((2.0, 3.0), (0.0, 3.0), (0.0, 0.0))
            value, g = gradient(term_product, a, b, x)
            @test value == a * b
            @test g[2:3] == (b, a)
        end
        for (a, b) in ((2.0, 3.0), (1.0, -1.0), (0.0, 0.0))
            value, g = gradient(term_sum, a, b, x)
            @test value == a + b
            @test g[2:3] == (1.0, 1.0)
        end
        for (f, value, expected) in (
            (poly_derivative, 4.0, 2.0),
            (term_constructor, 2.0, 1.0),
            (poly_constructor, 2.0, 1.0),
            (unary, -2.0, -1.0),
            (left_scalar, -2.0, 0.0),
            (integer_constant, 4.0, 1.0),
        )
            v, g = gradient(f, 2.0, x)
            @test v == value
            @test g[2] == expected
        end
        v, g = gradient(poly_scalar, 2.0, 3.0, x)
        @test v == 2.0
        @test collect(g[2:3]) ≈ [1.0, 0.0]
        for n in (0, 1, 3), a in (0.0, 2.0)
            v, g = gradient(poly_power, a, x, n)
            @test v == n * a
            @test g[2] == n
            v, g = gradient(term_power, a, x, n)
            @test v == a^n
            @test g[2] == (iszero(n) ? 0.0 : n * a^(n - 1))
        end
    end
    @testset "Coefficient arrays and cancelled support" begin
        cache = MC.prepare_gradient_cache(poly_coefficient, ones(3), x)
        for c in ([1.0, 2.0, 3.0], zeros(3), [1.0, -1.0, 3.0])
            before = copy(c)
            v, g = MC.value_and_gradient!!(cache, poly_coefficient, c, x)
            @test v == c[2]
            @test g[2] == [0.0, 1.0, 0.0]
            @test c == before
        end
        for c in ([2.0, 3.0, 4.0, 5.0], [1.0, 1.0, 1.0, -1.0], zeros(4))
            v, g = gradient(poly_product, c, x)
            @test v == c[1] * c[4] + c[2] * c[3]
            @test g[2] == c[[4, 3, 2, 1]]
        end
    end
    @testset "Coefficient views and Float32" begin
        for T in (Float32, Float64)
            c = T[1, 2, 3, 4]
            v, g = gradient(view_coefficient, c, x)
            @test v === T(2)
            @test g[2] == T[0, 1, 0, 0]
        end
    end
    @testset "Complex coefficient pairing" begin
        c = ComplexF64[2+im, 3-im, -1+2im, 4+3im]
        v, g = gradient(complex_product, c, x)
        @test v == complex_product(c, x)
        @test g[2] == conj.(c[[4, 3, 2, 1]])
        # Real and imaginary perturbations independently check conjugation.
        for i in eachindex(c), direction in (1.0, 1.0im)
            h = zeros(ComplexF64, length(c))
            h[i] = 1e-5 * direction
            fd = (complex_product(c + h, x) - complex_product(c - h, x)) / 2e-5
            @test real(conj(direction) * g[2][i]) ≈ fd rtol = 1e-8 atol = 1e-8
        end
        _, g = gradient(complex_mixed, 2.0, x)
        @test g[2] == 0.0
        # The scalar pullback uses SA's coefficient dot, with conjugation
        # on the polynomial's coefficients rather than the cotangent.
        p = (1.0 + 2im) * x + 3.0 * y
        v, g = gradient(scaled_coefficient, 3.0 + 4im, p, x)
        @test v == -5.0
        @test g[2] == 1.0 - 2im
        @test g[3].polynomial == (3.0 - 4im) * x
        @test p == (1.0 + 2im) * x + 3.0 * y
    end
    @testset "Polynomial inputs and aliasing" begin
        p, q = MP.polynomial(2.0 * x), MP.polynomial(3.0 * y)
        v, g = gradient(missing_coefficient, p, x)
        @test iszero(v)
        @test g[2].polynomial == x^3
        v, g = gradient(product_different_parents, p, q, x, y)
        @test v == 6.0
        @test g[2].polynomial == 3x
        @test g[3].polynomial == 2y
        v, g = gradient(repeated_argument, p, x)
        @test v == 4.0
        @test g[2].polynomial == 4x + x^2
        v, g = gradient(product_same_argument, p, p, x)
        @test v == 4.0
        # Mooncake prepares a separate gradient for each positional input.
        @test g[2].polynomial == 2x
        @test g[3].polynomial == 2x
        # Aliasing within a single input must share one tangent.
        tangents = MC.zero_tangent((p, p))
        @test tangents[1] === tangents[2]
        v, g = gradient(identity_polynomial, p, x)
        @test v == 2.0
        @test g[2].polynomial == x
        @test p == 2x
        @test q == 3y
    end
    @testset "Mooncake tangent interface" begin
        rng = MersenneTwister(0)
        for p in (
            MP.polynomial(2.0x),
            MP.polynomial((2.0 + im) * x),
            zero(MP.polynomial(2.0x)),
        )
            MC.TestUtils.test_tangent_interface(rng, p)
            MC.TestUtils.test_tangent_splitting(rng, p; test_opt_flag = false)
            MC.TestUtils.test_tangent_interface(rng, (p, p))
        end
    end
end

@testset "Mooncake" begin
    @test Ext !== nothing
    @static if isdefined(Mod, Symbol("@ncpolyvar"))
        Mod.@ncpolyvar nx ny
        @test_throws ArgumentError CRC.rrule(*, 1.0nx, 1.0ny)
        @test_throws ArgumentError CRC.rrule(^, 1.0nx + 1.0ny, 2)
    end
    Mod.@polyvar x y
    test_backend(x, y)
end

end

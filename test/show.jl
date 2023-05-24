struct CustomLaTeXPrint end

Base.:-(::CustomLaTeXPrint) = CustomLaTeXPrint()
Base.iszero(::CustomLaTeXPrint) = false
Base.show(io::IO, ::MIME"text/latex", ::CustomLaTeXPrint) = print(io, " \$\$ \\[\\(a_a \\) \\]\t  \$\$")

@testset "Show" begin
    Mod.@polyvar x y z
    @test sprint(show, (x * y^2 + x + 1 + y)) == "1 + y + x + xy²"
    @test sprint(show, (x + 1 + y) / x^2) == "(1 + y + x) / (x²)"
    @test sprint(show, (x - y - x + y) / (x^2 - x)) == "(0) / (-x + x²)"
    cp = CustomPoly(1 + x)
    @test sprint(show, cp) == "$(typeof(cp))(1 + x)"
    # Test taken from TypedPolynomials
    @test sprint(show, x) == "x"
    @test sprint(show, x^0) == "1"
    @test sprint(show, x^1) == "x"
    @test sprint(show, x^2) == "x²"
    @test sprint(show, 1x^2) == "x²"
    @test sprint(show, 5x) == "5x"
    @test sprint(show, x * y) == "xy"
    @test sprint(show, y * x) == "xy"
    @test sprint(show, 5 + y + x) == "5 + y + x"
    @test sprint(show, y + 5 + x) == "5 + y + x"
    @test sprint(show, x + x^2) == "x + x²"
    @test sprint(show, x^2 + x) == "x + x²"
    @test sprint(show, x^2 - 3.0x) == "-3.0x + x²"
    @test sprint(show, -x^3) == "-x³"
    @test sprint(show, -x^3 + y^2 - x) == "-x + y² - x³"
    @test sprint(show, -2.0x^2) == "-2.0x²"
    @test sprint(show, (1.0 + 3.1im) * x * z) == "(1.0 + 3.1im)xz"
    @test sprint(show, -(1.0 + 3.1im) * z * x) == "(-1.0 - 3.1im)xz"
    @test sprint(show, x^2 + (1.0 + 3.1im) * x) == "(1.0 + 3.1im)x + x²"
    @test sprint(show, x^2 - (1.0 + 3.1im) * x) == "(-1.0 - 3.1im)x + x²"

    @test sprint(show, 1e-8 + 3.2e7 * x + 4.5 * x^2) ==
          "1.0e-8 + 3.2e7x + 4.5x²"
    @test sprint(show, "text/plain", 1e-8 + 3.2e7 * x + 4.5 * x^2) ==
          "1.0e-8 + 3.2e7x + 4.5x²"
    @test sprint(show, "text/latex", 1e-8 + 3.2e7 * x + 4.5 * x^2) ==
          "\$\$ 1.0 \\cdot 10^{-8} + 3.2 \\cdot 10^{7}x + 4.5x^{2} \$\$"

    Mod.@polyvar x[0:9]
    @test sprint(show, sum(i * x[i]^i for i in 1:10)) ==
          "x₀ + 2x₁² + 3x₂³ + 4x₃⁴ + 5x₄⁵ + 6x₅⁶ + 7x₆⁷ + 8x₇⁸ + 9x₈⁹ + 10x₉¹⁰"
    @test sprint(show, "text/latex", sum(i * x[i]^i for i in 1:10)) ==
          "\$\$ x_{0} + 2x_{1}^{2} + 3x_{2}^{3} + 4x_{3}^{4} + 5x_{4}^{5} + 6x_{5}^{6} + 7x_{6}^{7} + 8x_{7}^{8} + 9x_{8}^{9} + 10x_{9}^{10} \$\$"
    @test sprint(show, "text/latex", (x[2] + 1) / x[3]^2) ==
          "\$\$ \\frac{1 + x_{1}}{x_{2}^{2}} \$\$"

    Mod.@polyvar x[1:11]
    @test sprint(show, "text/latex", x[10]) == "\$\$ x_{10} \$\$"
    @test sprint(show, x[10]) == "x₁₀"

    @test sprint(print, 2x[1]^2 + 3x[3] + 1 + x[4]) ==
          "1 + x[4] + 3*x[3] + 2*x[1]^2"

    a = CustomLaTeXPrint()
    @test sprint((io, x) -> show(io, "text/latex", x), term(a, x[1]^2)) ==
          "\$\$ (a_a)x_{1}^{2} \$\$"
    @test sprint(
        (io, x) -> show(io, "text/latex", x),
        polynomial([a, a], [x[1]^2, x[2]]),
    ) == "\$\$ (a_a)x_{2} + (a_a)x_{1}^{2} \$\$"
end

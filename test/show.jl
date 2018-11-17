@testset "Show" begin
    Mod.@polyvar x y z
    @test sprint(show, (x*y^2 + x + 1 + y)) == "xy² + x + y + 1"
    @test sprint(show, (x + 1 + y) / x^2) == "(x + y + 1) / (x²)"
    @test sprint(show, (x - y - x + y) / (x^2 - x)) == "(0) / (x² - x)"
    @test sprint(show, CustomPoly(1 + x)) == "CustomPoly{Int64,DynamicPolynomials.Polynomial{true,Int64}}(x + 1)"
    # Test taken from TypedPolynomials
    @test sprint(show, x) == "x"
    @test sprint(show, x^0) == "1"
    @test sprint(show, x^1) == "x"
    @test sprint(show, x^2) == "x²"
    @test sprint(show, 1x^2) == "x²"
    @test sprint(show, 5x) == "5x"
    @test sprint(show, x * y) == "xy"
    @test sprint(show, y * x) == "xy"
    @test sprint(show, 5 + y + x) == "x + y + 5"
    @test sprint(show, y + 5 + x) == "x + y + 5"
    @test sprint(show, x + x^2) == "x² + x"
    @test sprint(show, x^2 + x) == "x² + x"
    @test sprint(show, x^2 - 3.0x) == "x² - 3.0x"
    @test sprint(show, -x^3) == "-x³"
    @test sprint(show, -x^3 + y^2 - x) == "-x³ + y² - x"
    @test sprint(show, -2.0x^2) == "-2.0x²"
    @test sprint(show, (1.0 + 3.1im) * x*z) == "(1.0 + 3.1im)xz"
    @test sprint(show, -(1.0 + 3.1im) * z*x) == "(-1.0 - 3.1im)xz"
    @test sprint(show, x^2 + (1.0 + 3.1im) * x) == "x² + (1.0 + 3.1im)x"
    @test sprint(show, x^2 - (1.0 + 3.1im) * x) == "x² + (-1.0 - 3.1im)x"
    @test sprint(show, [1.0, 2.0] * x) == "([1.0, 2.0])x"

    Mod.@polyvar x[0:9]
    @test sprint(show, sum(i*x[i]^i for i=1:10)) == "10x₉¹⁰ + 9x₈⁹ + 8x₇⁸ + 7x₆⁷ + 6x₅⁶ + 5x₄⁵ + 4x₃⁴ + 3x₂³ + 2x₁² + x₀"
    @test sprint(show, "text/latex", sum(i*x[i]^i for i=1:10)) == "\$\$ 10x_{9}^{10} + 9x_{8}^{9} + 8x_{7}^{8} + 7x_{6}^{7} + 6x_{5}^{6} + 5x_{4}^{5} + 4x_{3}^{4} + 3x_{2}^{3} + 2x_{1}^{2} + x_{0} \$\$"
    @test sprint(show, "text/latex", (x[2] + 1) / x[3]^2) == "\$\$ \\frac{x_{1} + 1}{x_{2}^{2}} \$\$"

    Mod.@polyvar x[1:11]
    @test sprint(show, "text/latex", x[10]) == "\$\$ x_{10} \$\$"
    @test sprint(show, x[10]) == "x₁₀"

    @test sprint(print, 2x[1]^2+3x[3]+1+x[4]) == "2*x[1]^2 + 3*x[3] + x[4] + 1"
end

@testset "Show" begin
    Mod.@polyvar x y z
    @test sprint(show, (x*y^2 + x + 1 + y)) == "xy^2 + x + y + 1"
    @test sprint(show, (x + 1 + y) / x^2) == "(x + y + 1) / (x^2)"
    @test sprint(show, (x - y - x + y) / (x^2 - x)) == "(0) / (x^2 - x)"
    # Test taken from TypedPolynomials
    @test sprint(show, x) == "x"
    @test sprint(show, x^0) == "1"
    @test sprint(show, x^1) == "x"
    @test sprint(show, x^2) == "x^2"
    @test sprint(show, 1x^2) == "x^2"
    @test sprint(show, 5x) == "5x"
    @test sprint(show, x * y) == "xy"
    @test sprint(show, y * x) == "xy"
    @test sprint(show, 5 + y + x) == "x + y + 5"
    @test sprint(show, y + 5 + x) == "x + y + 5"
    @test sprint(show, x + x^2) == "x^2 + x"
    @test sprint(show, x^2 + x) == "x^2 + x"
end

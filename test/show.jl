@testset "Show" begin
    @polyvar x y
    @test sprint(show, SOSDecomposition([x+y, x-y])) == "(x + y)^2 + (x + -1y)^2"
    @test sprint(show, (x + 1 + y) / x^2) == "(x + y + 1) / (x^2)"
end

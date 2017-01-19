@testset "Promotion" begin
    @polyvar x y
    @inferred x*y+x
    @test typeof([x, x*y+x, x]) == Vector{Polynomial{true, Int}}
    @test typeof([1, x/y, x]) == Vector{RationalPoly{true, Int, Int}}
    @test typeof([(x^2-2x+2) x; x x^2]) == Matrix{Polynomial{true, Int}}
    @test typeof([2.0x, 3x]) == Vector{Term{true, Float64}}
    @inferred Any[x*y, x+y]
    @test typeof(Any[x*y, x+y]) == Vector{Any}
    @test typeof([x*y, x+y]) == Vector{Polynomial{true, Int}}
    @test typeof([2.0x, x/y, 1y]) == Vector{RationalPoly{true, Float64, Int}}
    @test typeof([2x+y, x/2.0y, x+1y]) == Vector{RationalPoly{true, Int, Float64}}

    X = [x, y]
    Y = [1 2; 3 4] * X
    @test Y[1] == x + 2y
    @test Y[2] == 3x + 4y
    Y = X' * [1 2; 3 4]
    @test Y[1] == x + 3y
    @test Y[2] == 2x + 4y
    @test dot(X, [1 2; 3 4] * X) == x^2 + 5x*y + 4y^2
end

@testset "Algebra" begin
    Mod.@polyvar x y
    @test 2 .- ((1 .+ (-x)) .* 4) ./ 2 == x.^2 .* (1 ./ x) .* 2
    @test dot(0, x^2 - 2*x^2) == dot((x^2 - x)', x^2 - x^2)
    @test 2 .* x .+ 2 == (x + 3) .+ (x .- 1)
    @test ((x + y) .- y) ./ y == x / y
    @test (-2)*x == -(2*x)
    @test x * x == x^2
    @test 4x == 2.0(2x)
    @test x * (2x) == 2x^2
    @inferred 2 * (2x)
    @inferred 2.0 * (2x)
    #@test eltype(2.0 * (2x)) == Float64
    @inferred 1.5 * (2x)
    @test 1.5 * (2x) isa AbstractTerm{Float64}
    @test 1.5 * (2x) == 3x

    @testset "Inference" begin
        @inferred x^2
        @inferred x^2-2x
        @inferred x+2x
        @inferred -2x+2
        @inferred x^2-2x+2
        @inferred x^2-2x
        @inferred 2x+2
        @inferred (1+x)^0
        @inferred (1+x)^1
        @inferred (1+x)^2
        @inferred (1+x)^3
        @inferred (x^2)^0
        @inferred (x^2)^1
        @inferred (x^2)^2
        @inferred (x^2)^3
        @inferred (2x) * (3x)
        @inferred (2x)^0
        @inferred (2x)^1
        @inferred (2x)^2
        @inferred (2x)^3
    end

    @test iszero((x+x-2*x) * (x * (x^2 + y^2)))
    @test iszero((0*x) * (x*y * (x^2 + y^2)))
end

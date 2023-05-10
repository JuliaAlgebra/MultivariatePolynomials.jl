@testset "Algebra" begin
    Mod.@polyvar x y
    @test 2 .- ((1 .+ (-x)) .* 4) ./ 2 == x .^ 2 .* (1 ./ x) .* 2
    @test dot(0, x^2 - 2 * x^2) == dot((x^2 - x)', x^2 - x^2)
    @test dot(x + 1, 2) == 2x + 2
    @test 2 .* x .+ 2 == (x + 3) .+ (x .- 1)
    @test ((x + y) .- y) ./ y == x / y
    @test (-2) * x == -(2 * x)
    @test x * x == x^2
    @test 4x == 2.0(2x)
    @test x * (2x) == 2x^2
    @inferred 2 * (2x)
    @inferred 2.0 * (2x)
    #@test eltype(2.0 * (2x)) == Float64
    t = @inferred 1.5 * (2x)
    @test t isa AbstractTerm{Float64}
    @test t == 3x

    p = @inferred y^2 * (2.0x^2 + 3y)
    @test p == 2x^2 * y^2 + 3y^3
    @test p isa AbstractPolynomial{Float64}

    p = @inferred (2x^2 + 3.0y) * y^2
    @test p == 2x^2 * y^2 + 3y^3
    @test p isa AbstractPolynomial{Float64}

    p = @inferred (2.0y^2) * (2x^2 + 3y)
    @test p == 4x^2 * y^2 + 6y^3
    @test p isa AbstractPolynomial{Float64}

    p = @inferred (2x^2 + 3y) * (2.0y^2)
    @test p == 4x^2 * y^2 + 6y^3
    @test p isa AbstractPolynomial{Float64}

    p = @inferred (y^2 + 2.0) * (2x + y)
    @test p == 2x * y^2 + y^3 + 4x + 2y
    @test p isa AbstractPolynomial{Float64}

    p = @inferred CustomPoly(x + 2y) * CustomTerms(x^2 + 1.0)
    @test p == x^3 + 2x^2 * y + x + 2y
    @test p isa AbstractPolynomial{Float64}

    p = @inferred CustomPoly(x + 2y) + 1
    @test p == x + 2y + 1
    @test p isa AbstractPolynomial{Int}

    p = @inferred CustomPoly(x + 2y) + 1.0
    @test p == x + 2y + 1
    @test p isa AbstractPolynomial{Float64}

    @testset "Inference" begin
        @inferred x^2
        @inferred x^2 - 2x
        @inferred x + 2x
        @inferred -2x + 2
        @inferred x^2 - 2x + 2
        @inferred x^2 - 2x
        @inferred 2x + 2
        @inferred (1 + x)^0
        @inferred (1 + x)^1
        @inferred (1 + x)^2
        @inferred (1 + x)^3
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

    @test iszero((x + x - 2 * x) * (x * (x^2 + y^2)))
    @test iszero((0 * x) * (x * y * (x^2 + y^2)))

    @testset "Scalar - Array" begin
        Mod.@polyvar x y
        @test x + [x^2+y y; y x*y] == [x+x^2+y x+y; x+y x+x*y]
        @test [x^2+y y; y x*y] + x == [x^2+x+y y+x; y+x x*y+x]

        @test x * [1 + y y] == [x * (1 + y) x * y]
        @test [1 + y y] * x == [(1 + y) * x y * x]
        @test [1 + y y] / x == [(1 + y) / x y / x]
    end
end

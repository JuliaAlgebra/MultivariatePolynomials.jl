@testset "Promotion" begin
    @polyvar x y
    @inferred x*y+x
    @test [x, x*y+x, x] isa Vector{<:AbstractPolynomial{Int}}
    @test [1, x/y, x] isa Vector{RationalPoly{<:AbstractTerm{Int}, <:AbstractVariable}}
    @test [(x^2-2x+2) x; x x^2] isa Matrix{<:AbstractPolynomial{Int}}
    @test [2.0x, 3x] == Vector{<:AbstractTerm{Float64}}
    @inferred Any[x*y, x+y]
    @test Any[x*y, x+y] isa Vector{Any}
    @test [x*y, x+y] == Vector{AbstractPolynomial{Int}}
    @test [2x*y, x+y] == Vector{AbstractPolynomial{Int}}
    @test [2.0x, x/y, 1y] == Vector{RationalPoly{<:AbstractTerm{Float64}, <:AbstractVariable}}
    @test [2x+y, x/2.0y, x+1y] == Vector{RationalPoly{<:AbstractPolynomial{Int}, <:AbstractTerm{Float64}}}

    X = [x, y]
    Y = [1 2; 3 4] * X
    @test Y[1] == x + 2y
    @test Y[2] == 3x + 4y
    Y = X' * [1 2; 3 4]
    @test Y[1] == x + 3y
    @test Y[2] == 2x + 4y
    @test dot(X, [1 2; 3 4] * X) == x^2 + 5x*y + 4y^2
end

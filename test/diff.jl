@testset "Differentiation" begin
    Mod.@polyvar x y
    @test differentiate(true*x+true*x^2, y) == 0
    @inferred differentiate(true*x+true*x^2, y)
    @test differentiate(x*y + 3y^2 , [x, y]) == [y, x+6y]
    @inferred differentiate(x*y + 3y^2 , x)
    @inferred differentiate(x*y + 3y^2 , [x, y])
    @test differentiate(1 / x , [x, y]) == [-1/x^2, 0]
    @test differentiate((x - y) / (x * y) , [x, y]) == [y^2 / (x * y)^2, -x^2 / (x * y)^2]
    p = x^2-y+x*y
    @test differentiate(p, x, 2) == 2
    @test differentiate(p, x, 0) === p
    @test_throws DomainError differentiate(p, x, -1)
    @test differentiate(x, x, 1) == 1
    #@inferred differentiate(x, x, 0) # FIXME failing at the moment
    @inferred differentiate(x, x, 1)
    @test differentiate(4x*y^3, y) == 12x*y^2
    @inferred differentiate(4x*y^3, y)
    @test differentiate(differentiate(x*y^4, y), y) == 12x*y^2
    @inferred differentiate(differentiate(x*y^4, y), y)
    @test differentiate(x*y^4, y, 3) == 24x*y
    @inferred differentiate(x*y^4, y, 3)
    @test differentiate(x*y^4, y, 0) == x*y^4
    @test differentiate(2x^2, x, 2) == 4
    @inferred differentiate(2x^2, x, 2)
    @inferred differentiate(2x^2 + x*y + y^2, [x, y])
    p = differentiate(2x^2 + 3x*y + y^2, [x, y], 2)
    @test isa(p, Matrix{<:AbstractPolynomial{Int}})
    @test p == [4 3; 3 2]
    p = differentiate(2x^2 + 3x*y^2 + 4y^3 + 2.0, (x, y), 2)
    @test isa(p, Matrix{<:AbstractPolynomial{Float64}})
    @test p == [4.0 6.0y; 6.0y 6.0x+24.0y]
end

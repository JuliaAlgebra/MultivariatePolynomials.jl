@testset "Norm of polynomial" begin
    @polyvar x y
    @test norm(x, 2)  == 1.0
    @test norm(2x, 2) == 2.0
    @test norm(2x+1, 2) ≈ sqrt(5)
    @test norm(2x*y+x^2+3) ≈ sqrt(14)
end

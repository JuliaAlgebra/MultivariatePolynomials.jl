@testset "Norm of polynomial" begin
    Mod.@polyvar x y
    @test norm(x) == norm(x, 2)  == 1.0
    @test norm(2x) == norm(2x, 2) == 2.0
    @test norm(2x+1) == norm(2x+1, 2) ≈ sqrt(5)
    @test norm(2x*y+x^2+3) ≈ norm(2x*y+x^2+3,2) ≈ sqrt(14)
end

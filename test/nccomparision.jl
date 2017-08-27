@testset "Comparison of monomials" begin
    @ncpolyvar x y
    @test x*y != y*x
end

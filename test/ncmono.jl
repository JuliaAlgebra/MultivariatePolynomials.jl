@testset "Non-commutative monomial" begin
    @ncpolyvar x y
    @test (x*x*y^2*x).z == [2, 2, 1]
end

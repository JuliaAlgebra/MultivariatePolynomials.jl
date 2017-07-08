@testset "Measure" begin
    @polyvar x y
    @test_throws ErrorException Measure([1, 2], [x, x*y, y])
    @test_throws ErrorException Measure([1, 2, 3, 4], [x, x*y, y])
    m = Measure([1, 0, 2, 3], [x^2*y^2, y*x^2, x*y*x^2, x*y^2])
    @test m.a == [2, 1, 0, 3]
end

@testset "Measure" begin
    @polyvar x y
    @test_throws ErrorException Measure([1, 2], [x, x*y, y])
    @test_throws ErrorException Measure([1, 2, 3, 4], MonomialVector([x, x*y, y]))
    m = Measure([1, 0, 2, 3], [x, x*y, y^4, x^2*y])
    @test m.a == [2, 3, 1]
end

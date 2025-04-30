@testset "Det" begin
    Mod.@polyvar x y

    @test det([x 1; y 1]) == x - y
    @test det([x 1; x 1]) == 0
    @test det([x 1 1 1; x y 1 2; 0 0 0 1; x 0 y 0]) == -x * y^2 + 2 * x * y - x
end

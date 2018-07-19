@testset "Det" begin
    Mod.@polyvar x y

    @test det([x 1; y 1]) == x - y
    @test det([x 1; x 1]) == 0
end

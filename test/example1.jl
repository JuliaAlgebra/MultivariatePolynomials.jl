@testset "Example 1" begin
    Mod.@polyvar x y
    p = 2x + 3.0x * y^2 + y
    @test differentiate(p, x) == 2 + 3y^2
    @test differentiate.(p, (x, y)) == (2 + 3y^2, 6x * y + 1)
    @test p((x, y) => (y, x)) == 2y + 3y * x^2 + x
    @test subs(p, y => x^2) == 2x + 3x^5 + x^2
    @test p(x => 1, y => 2) == 2 + 3 * 4 + 2
end

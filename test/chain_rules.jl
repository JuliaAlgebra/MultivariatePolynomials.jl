using Test
using ChainRulesCore

@testset "ChainRulesCore" begin
    Mod.@polyvar x y
    p = x + y
    q = x - y
    output, pullback = ChainRulesCore.rrule(+, p, q)
    @test output == 2x
    @test pullback(2) == (NoTangent(), 2, 2)
    @test pullback(x + 3) == (NoTangent(), x + 3, x + 3)
end

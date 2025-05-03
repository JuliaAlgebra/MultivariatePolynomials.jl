if !@isdefined(Mod)
    import DynamicPolynomials as Mod
    include("../utils.jl")
end

@testset "Comparison of monomials" begin
    Mod.@ncpolyvar x y
    @test x * y != y * x
    @test x^2 == x * x
end

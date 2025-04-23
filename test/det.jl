@testset "Det" begin
    Mod.@polyvar x y

    @testset "zero-one $T" for T in (Bool, Int, Float64, BigInt, BigFloat)
        z = zero(T)
        o = one(T)
        @test (@inferred det([x o; y o])) == x - y
        @test (@inferred det([x o; x o])) == 0
        @test (@inferred det([x+y y; z y])) == (x + y)*y
    end

    @test (@inferred det([x 1 1 1; x y 1 2; 0 0 0 1; x 0 y 0])) ==
          -x * y^2 + 2 * x * y - x
end

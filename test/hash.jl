@testset "Hashing" begin
    Mod.@polyvar x y z
    @test hash(x) != hash(y)
    @test hash(1x) == hash(1.0x)
    @test hash(1x + 3y) == hash(1.0x + 3.0y)
    @test hash(one(x)) == hash(x^0)
    @test hash(x * y) == hash(polynomial(x * y))
    @test hash(0.0x) == hash(0.0)
    @test hash(1.0x) == hash(x)
    @test hash(x - x) == hash(zero(x))
    @test hash(monomial_vector([z^3, z * x, y])[2:3]) ==
          hash(monomial_vector([z^3, z * x]))
    @test hash(monomial_vector([z^3, x, y])[2:2]) ==
          hash(monomial_vector([z^3, x])[1:1])
    @test hash(empty_monomial_vector(x)) == hash([])
    @test hash(1) == hash(one(x))
    @test hash(1) == hash(constant_monomial(x * y))
    @test hash(2) != hash(constant_monomial(x * y))
    @test hash(0.0) == hash(x - x)
end

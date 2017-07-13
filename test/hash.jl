@testset "Hashing" begin
    @polyvar x y z
    @test hash(x) != hash(y)
    @test hash(1x) == hash(1.0x)
    @test hash(1x+3y) == hash(1.0x+3.0y)
    @test hash(one(x)) == hash(x^0)
    @test hash(x*y) == hash(polynomial(x*y))
    @test hash(Term(1.0, Monomial(x))) == hash(x)
    @test hash(x-x) == hash(zero(x))
    @test hash(MonomialVector([z, y, x], [[3, 0, 0], [1,0,1]])) == hash(MonomialVector([z, x], [[3, 0], [1, 1]]))
    @test hash(Monomial([z, y, x], [3, 0, 0])) == hash(Monomial([z, x], [3, 0]))
    @test hash(1) == hash(one(x))
    @test hash(1) == hash(Monomial([x, y], [0, 0]))
    @test hash(2) != hash(Monomial([x, y], [0, 0]))
    @test hash(0.0) == hash(x-x)
end

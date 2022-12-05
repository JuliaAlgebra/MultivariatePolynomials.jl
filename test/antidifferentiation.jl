@testset "Antidifferentiation" begin
    Mod.@polyvar x y z
    @test antidifferentiate(3, y) == 3y
    @test antidifferentiate.([x, y], y) == [x * y, y^2 / 2]
    @test antidifferentiate(true * x + true * x^2, y) == x * y + x^2 * y
    @inferred antidifferentiate(true * x + true * x^2, y)
    @test antidifferentiate(x * y + 3y^2, [x, y]) == [x^2 * y / 2 + 3 * x * y^2, x * y^2 / 2 + y^3]
    @inferred antidifferentiate(x * y + 3y^2, x)
    @inferred antidifferentiate(x * y + 3y^2, [x, y])
    p = x^2 - y + x * y
    @test antidifferentiate(p, x, 2) == x^4 / 12 - x^2 * y / 2 + x^3 * y / 6
    @test_throws DomainError antidifferentiate(p, x, -1)
    @test antidifferentiate(x, x, 1) == x^2 / 2
    @inferred antidifferentiate(x, x, 1)
    @test antidifferentiate(4x * y^3, y) == x * y^4
    @inferred antidifferentiate(4x * y^3, y)
    @test antidifferentiate(antidifferentiate(x * y^4, y), y) == x * y^6 / 30
    @inferred antidifferentiate(antidifferentiate(x * y^4, y), y)
    @test antidifferentiate(7 * 6 * 5 * x * y^4, y, 1) == 7 * 6 * x * y^5
    @test antidifferentiate(7 * 6 * 5 * x * y^4, y, 2) == 7 * x * y^6
    @test antidifferentiate(7 * 6 * 5 * x * y^4, y, 3) == x * y^7
    @inferred antidifferentiate(x * y^4, y, 3)
    @test antidifferentiate(x * y^4, y, 0) == x * y^4
    @test antidifferentiate(2x^2, x, 2) == x^4 / 6
    @inferred antidifferentiate(2x^2, x, 2)
    @inferred antidifferentiate(2x^2 + x * y + y^2, [x, y])

    @testset "antidifferentiation with Val{}" begin
        @test @inferred(antidifferentiate(x, x, Val{0}())) == x
        @test @inferred(antidifferentiate(x, x, Val{1}())) == x^2 / 2
        @test @inferred(antidifferentiate(x^2, x, Val{1}())) == x^3 / 3
        @test @inferred(antidifferentiate(x^2, y, Val{1}())) == x^2 * y
        @test @inferred(antidifferentiate(2x^2 + 3y, x, Val{1}())) == 2 / 3 * x^3 + 3 * x * y
        @test @inferred(antidifferentiate(2x^2 + 3y, x, Val{2}())) == x^4 / 6 + 3 / 2 * x^2 * y
    end
end

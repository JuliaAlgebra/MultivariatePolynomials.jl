@testset "Term" begin
    Mod.@polyvar x
    @test Any(1x) == 1x
    @test one(1x) == one(1.0x) == 1
    @test zero(1x) == zero(1.0x) == 0
    @test nvariables(0.0x) == 1
    @test nvariables(1x) == 1
    #@inferred one(1x)
    @inferred zero(1x)
    #@inferred one(1.0x)
    @inferred zero(1.0x)

    @test monic(2.0x) isa AbstractTerm{Float64}
    @test monic(2.0x) == x
    @test monic(2x^2) isa AbstractTerm{Int}
    @test monic(2x^2) == x^2

    @test leadingterm(2x^2) == 2x^2
    @test nterms(2x^2) == 1
    @test terms(2x^2) == [2x^2]
    @test nterms(0*x) == 0
    @test terms(0*x) == typeof(0*x)[]
    @test nterms(0.0x) == 0
    @test terms(0.0x) == typeof(0.0x)[]

    @test nterms(polynomial(0.0x)) == 0
    @test nterms(convert(polynomialtype(0.0x), 0.0x)) == 0

    @test term(x) isa AbstractTerm
    @test term(x^2) == x^2
    @test term(1x^2) isa AbstractTerm
    @test term(1x) == x
    @test zeroterm(1x) == 0*x

    Mod.@polyvar y
    @test degree(2x^2, x) == 2
    @test degree(2x^2, y) == 0
    @test degree(2x^2, y) == 0

    @test_throws InexactError push!([1], 2x)
    @test_throws ErrorException push!([x^2], 2x)
end

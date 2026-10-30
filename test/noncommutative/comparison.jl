if !@isdefined(Mod)
    import DynamicPolynomials as Mod
    include("../utils.jl")
end

@testset "Comparison of monomials" begin
    Mod.@ncpolyvar x y
    @test x * y != y * x
    @test x^2 == x * x
end

@testset "Noncommutative comparison overrides" begin
    Mod.@ncpolyvar x[1:3]
    monos = monomials(x, 2)
    for a in x, b in x
        word = a * b
        @test count(isequal(word), monos) == 1
        @test count(m -> iszero(cmp(word, m)), monos) == 1
        @test word in Set(monos)
    end

    # The default implementations require commutative monomials. A frontend
    # must supply its own methods to compare and hash noncommutative words.
    word = x[2] * x[1]
    signature = Tuple{AbstractMonomial,AbstractMonomial}
    for op in (==, isequal, cmp)
        @test_throws AssertionError invoke(op, signature, word, word)
    end
    @test_throws AssertionError invoke(
        hash,
        Tuple{AbstractMonomial,UInt},
        word,
        UInt(0),
    )
    @test_throws AssertionError invoke(
        hash,
        Tuple{MP.BasisElement,UInt},
        word,
        UInt(0),
    )
end

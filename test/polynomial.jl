const MP = MultivariatePolynomials
@testset "Polynomial" begin
    Mod.@polyvar x

    @test terms(polynomial([1, x^2, x, 2x^2])) == [3x^2, x, 1]
    @test terms(polynomial([x, 3x^4, 2], MP.UniqState())) == [3x^4, x, 2]
    @test terms(polynomial([x^3, 2x^3, x^2, -2x^2, x^2, x, 2, -2], MP.SortedState())) == [3x^3, x]

    @test polynomial(1 + x) == 1 + x
    @test leadingterm(1 + x) == x
    @test leadingterm(x - x) == 0
    @test leadingmonomial(x - x) == 1
    @test leadingcoefficient(x - x) == 0
    @test one(1 + x) == one(1.0 + x) == 1
    @test zero(1 + x) == zero(1.0 + x) == 0
    @test 1 != 1 + x
    @test 2x == x + x^2 + x - x^2
    @test x + x^2 - x^2 - x != 2x
    @test x^2 + x != x^2 + x + 1
    #@inferred one(1 + x)
    @inferred zero(1 + x)
    #@inferred one(1.0 + x)
    @inferred zero(1.0 + x)

    @test polynomial(1:2, [x, x^2]) == x + 2x^2
    @test polynomial(1:2, monomials(x, 1:2)) == 2x + x^2

    @test terms(polynomial(1 + x + x^2 - x + x^2)) == [2x^2, 1]
    @test terms(CustomPoly(1 + x + x^2 - x + x^2)) == [2x^2, 1]

    @test (1.0 + x) * x == x^2 + x
    @test constantterm(1, x) * (1 - x) == 1 - x
    @test promote_type(typeof(1-x), typeof(x)) <: AbstractPolynomial{Int}
    @test x != 1 - x

    @test term(x + x^2 - x) isa AbstractTerm
    @test term(x + x^2 - x) == x^2
    @test term(x - x) isa AbstractTerm
    @test iszero(term(x - x))
    @test_throws ErrorException term(x + x^2)

    Mod.@polyvar y

    p = 3x^2*y^4 + 2x
    @test terms(p)[end] == 2x
    typetests(p)
    typetests([p, x + y])
    @test (@inferred polynomial(p)) isa AbstractPolynomial{Int}
    @test (@inferred polynomial(p, Float64)) isa AbstractPolynomial{Float64}

    @test coefficient(2x + 4y^2 + 3, y^2) == 4
    @test coefficient(2x + 4y^2 + 3, x^2) == 0

    @test (@inferred 2x^2*y + 0.0x*y) == 2x^2*y
    @test (@inferred 0.0x^2*y + 3x*y) == 3x*y

    @test iszero(((x + x) - 2x) * (x * (x ^ 2 + y ^ 2)))

    @test Tuple(variables([x + 1, y^2])) == (x, y)
    @test Tuple(variables([y^2, x + 1])) == (x, y)

    @test maxdegree(x^2 - x^2) == 0
    @test maxdegree(x^2 - x^2, x) == 0
    @test maxdegree(x^2 - x^2, y) == 0
    @test mindegree(x^2 - x^2) == 0
    @test mindegree(x^2 - x^2, x) == 0
    @test mindegree(x^2 - x^2, y) == 0
    @test extdegree(x^2 - x^2) == (0, 0)
    @test extdegree(x^2 - x^2, x) == (0, 0)
    @test extdegree(x^2 - x^2, y) == (0, 0)
    @test maxdegree(x*y + 2 + x^2*y + x + y) == 3
    @test maxdegree(x*y + 2 + x^2*y + x + y, x) == 2
    @test maxdegree(x*y + 2 + x^2*y + x + y, y) == 1
    @test mindegree(x*y + 2 + x^2*y + x + y) == 0
    @test mindegree(x*y + 2 + x^2*y + x + y, x) == 0
    @test mindegree(x*y + 2 + x^2*y + x + y, y) == 0
    @test extdegree(x*y + 2 + x^2*y + x + y) == (0, 3)
    @test extdegree(x*y + 2 + x^2*y + x + y, x) == (0, 2)
    @test extdegree(x*y + 2 + x^2*y + x + y, y) == (0, 1)
    @test extdegree(x*y + x^2*y, x) == (1, 2)
    @test extdegree(x*y + x^2*y, y) == (1, 1)
    @test leadingterm(x*y + 2 + x^2*y + x + y) == x^2*y
    @test nvariables(x + y - x) == 2
    @test nvariables(x + x^2) == 1

    @test coefficients(x*y + 2 + 3x^2*y + 4x + 6y, [x, x*y^2, x*y, x^2*y, y, x^3]) == [4, 0, 1, 3, 6, 0]

    # Doc examples
    @test coefficients(4x^2*y + x*y + 2x) == [4, 1, 2]
    @test coefficients(4x^2*y + x*y + 2x + 3, [x, 1, x*y, y]) == [2, 3, 1, 0]

    p = polynomial([4, 9], [x, x*x])
    @test coefficients(p) == [9, 4]
    @test monomials(p)[1] == x^2
    @test monomials(p)[2] == x
    @test p == dot([4, 9], [x, x*x])

    @inferred polynomial(i -> float(i), [x, x*x])
    @inferred polynomial(i -> 3 - float(i), monovec([x*x, x]))
    for p in (polynomial(i -> float(i), [x, x*x]),
              polynomial(i -> 3 - float(i), monovec([x*x, x])))
        @test coefficients(p) == [2.0, 1.0]
        @test monomials(p) == monovec([x^2, x])
    end

    @test (x + y)' == x + y
    @test transpose(x + y) == x + y
    @test transpose([1 2; 3 4] * x) == [1 3; 2 4] * x

    @test removemonomials(4x^2*y + x*y + 2x, [x*y]) == 4x^2*y + 2x

    @test_throws InexactError push!([1], x+1)

    @test polynomial([1 2; 3 4], [x^2, y]) == x^4 + 5x^2*y + 4y^2
    @test polynomial([1 2; 3 4], [x^2, y], Float64) isa AbstractPolynomial{Float64}
    @test polynomial([1 2; 3 4], [y, x^2]) == y^2 + 5x^2*y + 4x^4
    @test polynomial([1 2; 3 4], [y, x^2], Float64) isa AbstractPolynomial{Float64}
    @test polynomial([1 2; 3 4], monovec([y, x^2])) == x^4 + 5x^2*y + 4y^2
    @test polynomial([1 2; 3 4], monovec([y, x^2]), Float64) isa AbstractPolynomial{Float64}

    @test (@inferred round(2.6x + 1.001x^2)) == 3x + 1x^2
    @test (@inferred round(3.1x*y)) == 3x*y
    @test (@inferred round(2.613x + 1.1051x^2, digits=2)) ≈ 2.61x + 1.11x^2
    @test (@inferred round(3.145x*y, digits=1)) ≈ 3.1x*y

    @testset "Graded Lex Order" begin
        Mod.@polyvar x y z
        p = 3*y^2 + 2*y*x
        @test coefficients(p) == [2, 3]
        @test monomials(p) == monovec([x*y, y^2])
        # Examples from p. 59 of the 4th edition of "Ideals, Varieties, and Algorithms" of Cox, Little and O'Shea
        f = 4*x*y^2*z + 4*z^2 - 5*x^3 + 7*x^2*z^2
        @test coefficients(f) == [7, 4, -5, 4]
        @test monomials(f) == monovec([x^2*z^2, x*y^2*z, x^3, z^2])
    end
end

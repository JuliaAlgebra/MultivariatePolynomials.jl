@testset "PolyVar" begin
    @polyvar x
    @test zero(x) == 0
    @test typeof(zero(x)) == VecPolynomial{Int}
    @inferred zero(x)
    @test one(x) == 1
    @test typeof(one(x)) == VecPolynomial{Int}
    @inferred one(x)
end
@testset "Monomial" begin
    @polyvar x
    @test_throws ArgumentError Monomial([x], [1,0])
    @test zero(x^2) == 0
    @test typeof(zero(x^2)) == VecPolynomial{Int}
    @inferred zero(x^2)
    @test one(x^2) == 1
    @test typeof(one(x^2)) == VecPolynomial{Int}
    @inferred one(x^2)
end
@testset "Degree to Monomial Vector" begin
    @polyvar x y z
    @test monomials([x,y], 2, z->z[1] <= 1) == [x*y,y^2]
    @test MonomialVector([x,y], 1:3) == MonomialVector([x^3,x^2*y,x*y^2,y^3,x^2,x*y,y^2,x,y])
    X = MonomialVector([x,y,z], 5)
    @test issorted(X, rev=true)
    @test X == MonomialVector([x^5, x^4*y, x^4*z, x^3*y^2, x^3*y*z, x^3*z^2, x^2*y^3, x^2*y^2*z, x^2*y*z^2, x^2*z^3, x*y^4, x*y^3*z, x*y^2*z^2, x*y*z^3, x*z^4, y^5, y^4*z, y^3*z^2, y^2*z^3, y*z^4, z^5])
end
@testset "Vector to Monomial Vector" begin
    @polyvar x
    @test MonomialVector([1]) == MonomialVector(PolyVar[], [Int[]])
    @test MonomialVector([x]) == MonomialVector([x], [[1]])
end
module newmodule
    using Base.Test
    import MultivariatePolynomials
    @testset "Polyvar macro hygiene" begin
        # Verify that the @polyvar macro works when the package has been activated
        # with `import` instead of `using`.
        MultivariatePolynomials.@polyvar x y
        @test isa(x, MultivariatePolynomials.PolyVar)
        @test isa(y, MultivariatePolynomials.PolyVar)
    end
end

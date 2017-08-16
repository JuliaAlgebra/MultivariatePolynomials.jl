import MultivariatePolynomials: similarvariable, @similarvariable

@testset "Create similar variable" begin
    Mod.@polyvar x y
    f = x^2 + y

    z = similarvariable(f, Val{:z})
    @test_broken z isa MultivariatePolynomials.AbstractVariable

    @inferred similarvariable(f, Val{:z})

    w = similarvariable(f, :w)
    @test_broken w isa MultivariatePolynomials

    @similarvariable f o
    @test_broken o isa MultivariatePolynomials.AbstractVariable
end

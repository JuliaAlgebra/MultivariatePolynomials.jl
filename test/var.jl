import MultivariatePolynomials: AbstractVariable, similarvariable, @similarvariable

@testset "Create similar variable" begin
    Mod.@polyvar x y
    f = x^2 + y

    z = similarvariable(f, Val{:z})
    @test z isa AbstractVariable

    @inferred similarvariable(f, Val{:z})

    w = similarvariable(f, :w)
    @test w isa AbstractVariable

    @similarvariable f o
    @test o isa AbstractVariable

    m = @similarvariable f u
    @test m isa AbstractVariable
end

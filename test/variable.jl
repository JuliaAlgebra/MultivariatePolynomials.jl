import MultivariatePolynomials: AbstractVariable, similarvariable, @similarvariable

@testset "Variable" begin
    @testset "polyvar macro index set" begin
        Mod.@polyvar x y z
        Mod.@polyvar x[1:3] y z[1:2]
        @test length(x) == 3
        @test length(z) == 2
        @test x[1] > x[2] > x[3] > y > z[1] > z[2]
    end
    @testset "PolyVar" begin
        Mod.@polyvar x
        alloc_test(() -> convert(typeof(x), x), 0)
        alloc_test(() -> convert(variable_union_type(x), x), 0)
        @test 1 != x
        @test x != 0
        @test copy(x) == x
        @test nvariables(x) == 1
        @test !isapproxzero(x)
        @test !iszero(x)
        @test zero(x) == 0
        @test iszero(zero(x))
        @test zero(x) isa AbstractPolynomial{Int}
        @inferred zero(x)

        typetests(x)
        @test (@inferred polynomial(x)) isa AbstractPolynomial{Int}
        @test (@inferred polynomial(x, Float64)) isa AbstractPolynomial{Float64}

        @test nterms(x) == 1
        @test @inferred(terms(x)) == [x]

        Mod.@polyvar y
        @test degree(x, x) == 1
        @test degree(x, y) == 0
        @test length(exponents(x)) == 1
        @test first(exponents(x)) == 1
        @test isconstant(x) == false

        @test divides(x, x) == true
        @test divides(x, y) == false

        @testset "Issue #82" begin
            for v in (x, typeof(x))
                for fun in (oneunit, one)
                    @test fun(v) == 1
                    @test 1 == fun(v)
                    #@test isone(fun(v)) # Enable when Julia v0.6 is dropped
                    @test fun(v) isa AbstractMonomial
                    @inferred fun(v)
                end
            end
        end

        @testset "Effective variables" begin
            @test [x] == @inferred effective_variables(x)
            @test [y] == @inferred effective_variables(y)
        end
    end
    @testset "Create similar variable" begin
        Mod.@polyvar x y
        f = x^2 + y

        z = similarvariable(f, Val{:z})
        @test z isa AbstractVariable

        z = similarvariable(typeof(f), Val{:z})
        @test z isa AbstractVariable

        @inferred similarvariable(f, Val{:z})
        @inferred similarvariable(typeof(f), Val{:z})

        w = similarvariable(f, :w)
        @test w isa AbstractVariable

        @similarvariable f o
        @test o isa AbstractVariable

        m = @similarvariable f u
        @test m isa AbstractVariable
    end
end

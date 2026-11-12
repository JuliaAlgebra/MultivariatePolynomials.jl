module TestBasisTypeDisplay

using Test
import MultivariatePolynomials as MP
import StarAlgebras as SA
include("SimplePolynomials.jl")

@testset "Compact basis type display" begin
    SimplePolynomials.@polyvar x y
    full = MP.FullBasis{MP.Monomial}([x, y])
    sub = SA.SubBasis(full, [[0, 0], [1, 0]])
    V = typeof(MP.variables(full))
    E = Vector{Int}
    for (basis, name) in ((full, "FullBasis"), (sub, "SubBasis"))
        expected = "MultivariatePolynomials.$name{$(MP.Monomial), $V, $E}"
        @test sprint(show, typeof(basis)) == expected
        @test sprint(show, Tuple{typeof(basis)}) == "Tuple{$expected}"
        @test sprint(show, Vector{typeof(basis)}) == "Vector{$expected}"
    end
    for (family, name) in
        ((MP.FullBasis, "FullBasis"), (MP.SubBasis, "SubBasis"))
        prefix = "MultivariatePolynomials.$name"
        @test sprint(show, family) == prefix
        @test sprint(show, family{MP.Monomial}) == "$prefix{$(MP.Monomial)}"
        @test sprint(show, family{MP.Monomial,V}) ==
              "$prefix{$(MP.Monomial), $V}"
        @test sprint(show, family{MP.Monomial,V,E}) ==
              "$prefix{$(MP.Monomial), $V, $E}"
        bounded = family{B} where {B<:MP.AbstractMonomialBasis}
        @test sprint(show, bounded) == "$prefix{<:$(MP.AbstractMonomialBasis)}"
        first_type = family{MP.Monomial,V,Vector{Int}}
        second_type = family{MP.Monomial,V,Vector{Int32}}
        displayed = sprint(show, Union{first_type,second_type})
        @test startswith(displayed, "Union{")
        @test occursin(sprint(show, first_type), displayed)
        @test occursin(sprint(show, second_type), displayed)
        lower_bound = family{MP.Monomial,V,TE} where {TE>:Vector{Int}}
        @test sprint(show, lower_bound) == sprint(
            (io, type) -> invoke(show, Tuple{IO,Type}, io, type),
            lower_bound,
        )
    end
    # Preserve variable names for implementations that encode them in a tuple.
    V = Tuple{Val{:x},Val{:y}}
    @test occursin(string(V), sprint(show, MP.FullBasis{MP.Monomial,V,E}))
    @test sprint(show, Union{}) == "Union{}"
    @test isempty(filter(m -> m.name == :show, Test.detect_unbound_args(MP)))
end

end # module

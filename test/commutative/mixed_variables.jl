module TestMixedVariables

using Test
import MultivariatePolynomials as MP

if !isdefined(Main, :Mod)
    include("../SimplePolynomials.jl")
end
const Mod = isdefined(Main, :Mod) ? Main.Mod : SimplePolynomials

# A separate instance supplies a distinct variable family without dependencies.
module OtherImplementation
include("../SimplePolynomials.jl")
end
const Other = OtherImplementation.SimplePolynomials

@testset "Incompatible variable families" begin
    Mod.@polyvar x y
    Other.@polyvar z
    A, B = MP.variable_union_type(x), MP.variable_union_type(z)
    @test A !== B
    err_xz = ArgumentError("Cannot combine variables of type $A and $B.")
    err_zx = ArgumentError("Cannot combine variables of type $B and $A.")
    for p in (x, x^2, 2x, x + 1), q in (z, z^2, 3z, z + 1)
        for op in (+, -, *)
            @test_throws err_xz op(p, q)
            @test_throws err_zx op(q, p)
        end
    end
    for op in (+, -, min, max)
        @test_throws err_xz MP.map_exponents(op, x, z)
        @test_throws err_zx MP.map_exponents(op, z, x)
    end
    vx, vz = MP.variables(x), MP.variables(z)
    # Empty variable collections still retain their backend.
    for (a, b) in ((vx, vz), (filter(_ -> false, vx), filter(_ -> false, vz)))
        @test_throws err_xz MP.promote_variables_with_maps(a, b)
        @test_throws err_zx MP.promote_variables_with_maps(b, a)
    end
    # TP variables may have different concrete types within the same family.
    @test MP.variable_union_type(x) === MP.variable_union_type(y)
    @test (@inferred x * y) == y * x
    @test (@inferred x + y) == y + x
    @test (@inferred x - y) == -(y - x)
end

end

module TestMixedVariables

using Test
import DynamicPolynomials as DP
import TypedPolynomials as TP
import MultivariatePolynomials as MP

@testset "Incompatible variable families (DynamicPolynomials#219)" begin
    DP.@polyvar x y
    DP.@ncpolyvar s
    TP.@polyvar z w
    for (a, b) in ((x, s), (x, z), (s, z))
        A, B = MP.variable_union_type(a), MP.variable_union_type(b)
        err_ab = ArgumentError("Cannot combine variables of type $A and $B.")
        err_ba = ArgumentError("Cannot combine variables of type $B and $A.")
        for p in (a, a^2, 2a, a + 1), q in (b, b^2, 3b, b + 1)
            for op in (+, -, *)
                @test_throws err_ab op(p, q)
                @test_throws err_ba op(q, p)
            end
        end
        for op in (+, -, min, max)
            @test_throws err_ab MP.map_exponents(op, a, b)
            @test_throws err_ba MP.map_exponents(op, b, a)
        end
        va, vb = MP.variables(a), MP.variables(b)
        # Empty variable collections still retain their backend.
        for (u, v) in
            ((va, vb), (filter(_ -> false, va), filter(_ -> false, vb)))
            @test_throws err_ab MP.promote_variables_with_maps(u, v)
            @test_throws err_ba MP.promote_variables_with_maps(v, u)
        end
    end
    # Different concrete TP variable types belong to the same family.
    @test typeof(z) !== typeof(w)
    @test MP.variable_union_type(z) === MP.variable_union_type(w)
    for (a, b) in ((x, y), (z, w))
        @test (@inferred a * b) == b * a
        @test (@inferred a + b) == b + a
        @test (@inferred a - b) == -(b - a)
    end
end

end

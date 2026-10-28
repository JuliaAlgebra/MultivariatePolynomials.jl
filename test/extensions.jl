module TestExtensions

using Test
if !isdefined(Main, :Mod)
    include("SimplePolynomials.jl")
end
const Mod = isdefined(Main, :Mod) ? Main.Mod : SimplePolynomials
import MultivariatePolynomials as MP
import JuMP
import MathOptInterface as MOI

function test_printing(var)
    x = MP.similar_variable(var, Val{:x})
    model = JuMP.Model()
    JuMP.@variable(model, a)
    for (p, plain, latex) in (
        (x, "x", "x"),
        (x^2, "x²", "x^{2}"),
        (2x^2, "2x²", "2x^{2}"),
        (1 + 2x^2, "1 + 2x²", "1 + 2x^{2}"),
        (a * x^2, "(a)x²", "(a)x^{2}"),
        (a * x^2 + 1, "(1) + (a)x²", "(1) + (a)x^{2}"),
    )
        @test JuMP.function_string(MIME"text/plain"(), p) == plain
        @test JuMP.function_string(MIME"text/latex"(), p) == latex
    end
    return
end

function test_without_indices(x)
    for op in (MOI.Utilities.map_indices, MOI.Utilities.substitute_variables)
        for m in (x, x^2), container in (m, (m, m), [m], [(m, m)], typeof(m)[])
            @test (@inferred op(error, container)) === container
        end
        for container in ((), Tuple{}[], (x, x^2), [(x, x^2)])
            @test (@inferred op(error, container)) === container
        end
    end
    return
end

@testset "JuMP and MOI extensions" begin
    Mod.@polyvar x
    test_printing(x)
    test_without_indices(x)
end

end

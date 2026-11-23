module TestInvalidations

using Test

# Invalidations can only be measured when the package is loaded for the first
# time so we run the measurement in a fresh process. The dependencies (and the
# triggers of the extensions) are loaded first so that only the invalidations
# caused by methods defined in MultivariatePolynomials (or its extensions) are
# counted.
const INVALIDATIONS_SCRIPT = """
using LinearAlgebra, MutableArithmetics, StarAlgebras, MathOptInterface, JuMP
using SnoopCompileCore
invs = @snoop_invalidations using MultivariatePolynomials
using SnoopCompile
dir = pkgdir(MultivariatePolynomials)
function belongs_to_mp(tree)
    if tree.method === nothing
        return false
    end
    return startswith(String(tree.method.file), dir)
end
@assert !belongs_to_mp((method = nothing,))
trees = filter(belongs_to_mp, invalidation_trees(invs))
show(stderr, MIME"text/plain"(), trees)
println(stderr)
for name in (:show, :other)
    selected = filter(trees) do tree
        return (tree.method.name == :show) == (name == :show)
    end
    println(sum(SnoopCompile.countchildren, selected; init = 0))
end
"""

@testset "Invalidations" begin
    project = dirname(Base.active_project())
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$project -e $INVALIDATIONS_SCRIPT`
    display_count, count = parse.(Int, split(readchomp(cmd)))
    # Adding type display methods invalidates callers of Base.show(::IO, ::Type)
    # compiled before MP was loaded. Track these separately to preserve the
    # tighter limit on arithmetic invalidations below.
    @test display_count <= (VERSION < v"1.11" ? 15_000 : 5_000)
    # The remaining invalidations are due to the equality between polynomials
    # and JuMP or MOI scalars defined by the extensions. They invalidate poorly
    # inferred call sites in MathOptInterface and JuMP:
    # `in(::Any, ::Set{MOI.VariableIndex})` in
    # `MOI.Utilities._build_copy_variables_with_set_cache`,
    # `==(::JuMP.NonlinearExpr, ::WeakRef)` and
    # `==(::Any, ::JuMP.NonlinearParameter)` in `JuMP.function_string`.
    # These will be fixed in JuMP and MOI, after which this should be 0.
    # The count depends on the Julia version, e.g., Julia v1.13.4 also records
    # `isequal(::JuMP.GenericNonlinearExpr, ::Any)` (81 -> 86) and Julia
    # v1.14.0-beta1 more callers of `==(::Any, ::JuMP.NonlinearParameter)`
    # (90 -> 95).
    @test count <= (VERSION >= v"1.14-" ? 95 : 86)
end

end

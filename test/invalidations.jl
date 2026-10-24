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
trees = filter(invalidation_trees(invs)) do tree
    return startswith(String(tree.method.file), dir)
end
show(stderr, MIME"text/plain"(), trees)
println(stderr)
println(sum(SnoopCompile.countchildren, trees; init = 0))
"""

@testset "Invalidations" begin
    project = dirname(Base.active_project())
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$project -e $INVALIDATIONS_SCRIPT`
    count = parse(Int, readchomp(cmd))
    # The remaining invalidations are due to the equality between polynomials
    # and JuMP or MOI scalars defined by the extensions. They invalidate poorly
    # inferred call sites in MathOptInterface and JuMP:
    # `in(::Any, ::Set{MOI.VariableIndex})` in
    # `MOI.Utilities._build_copy_variables_with_set_cache`,
    # `==(::JuMP.NonlinearExpr, ::WeakRef)` and
    # `==(::Any, ::JuMP.NonlinearParameter)` in `JuMP.function_string`.
    # These will be fixed in JuMP and MOI, after which this should be 0.
    @test count <= 81
end

end

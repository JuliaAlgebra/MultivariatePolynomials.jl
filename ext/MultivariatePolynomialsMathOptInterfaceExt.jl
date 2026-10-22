module MultivariatePolynomialsMathOptInterfaceExt

import MathOptInterface as MOI
import MultivariatePolynomials as MP

MP.@scalar_arithmetic MOI.AbstractScalarFunction

# Empty tuples are already handled by MOI, independently of this extension.
const ObjectOrTupleWithoutIndex = Union{
    MP.AbstractMonomialLike,
    Tuple{MP.AbstractMonomialLike,Vararg{MP.AbstractMonomialLike}},
}
const ObjectOrTupleOrArrayWithoutIndex =
    Union{ObjectOrTupleWithoutIndex,AbstractArray{<:ObjectOrTupleWithoutIndex}}

function MOI.Utilities.map_indices(
    ::F,
    x::ObjectOrTupleOrArrayWithoutIndex,
) where {F<:Function}
    return x
end

function MOI.Utilities.substitute_variables(
    ::F,
    x::ObjectOrTupleOrArrayWithoutIndex,
) where {F<:Function}
    return x
end

end

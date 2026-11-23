# These aliases live outside StarAlgebras, where the underlying types are defined.
# Explicit display methods also use the aliases inside other parametric types.
function _show_basis_type(io::IO, name, type)
    unwrapped = Base.unwrap_unionall(type)
    if type === Union{} || unwrapped isa Union
        return invoke(show, Tuple{IO,Type}, io, type)
    end
    element = unwrapped.parameters[1]
    parameters = Base.unwrap_unionall(element).parameters
    for parameter in parameters
        if parameter isa TypeVar && parameter.lb !== Union{}
            return invoke(show, Tuple{IO,Type}, io, type)
        end
    end
    last = length(parameters)
    while last > 0
        parameter = parameters[last]
        if !(parameter isa TypeVar && parameter.ub === Any)
            break
        end
        last -= 1
    end
    print(io, "MultivariatePolynomials.", name)
    if last > 0
        print(io, "{")
        for i in 1:last
            if i > 1
                print(io, ", ")
            end
            parameter = parameters[i]
            if parameter isa TypeVar
                print(io, "<:")
                show(io, parameter.ub)
            else
                show(io, parameter)
            end
        end
        print(io, "}")
    end
    return
end

# A single method avoids ambiguity at Union{}, which is below both families.
function Base.show(io::IO, type::Type{<:MonomialIndexedBasis})
    name = type <: FullBasis ? "FullBasis" : "SubBasis"
    return _show_basis_type(io, name, type)
end

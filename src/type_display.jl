# These aliases live outside StarAlgebras, where the underlying types are defined.
# Explicit display methods also use the aliases inside other parametric types.
function _show_basis_type(io::IO, name, type)
    if type isa Union
        return invoke(show, Tuple{IO,Type}, io, type)
    end
    element = Base.unwrap_unionall(type).parameters[1]
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

function Base.show(io::IO, type::Type{<:FullBasis})
    return _show_basis_type(io, "FullBasis", type)
end

function Base.show(io::IO, type::Type{<:SubBasis})
    return _show_basis_type(io, "SubBasis", type)
end

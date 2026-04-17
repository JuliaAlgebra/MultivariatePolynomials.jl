# From MultivariateBases — Variables struct

struct Variables{B,V}
    variables::V
end

Variables{B}(vars) where {B} = Variables{B,typeof(vars)}(vars)

function Base.one(v::Variables)
    return monomial(v.variables, constant_monomial_exponents(v))
end

function variable_index(v::Variables, var)
    return findfirst(isequal(var), v.variables)
end

function Base.:(==)(v::Variables{B}, w::Variables{B}) where {B}
    return v.variables === w.variables || v.variables == w.variables
end

monomial_type(::Type{Variables{B,V}}) where {B,V} = monomial_type(V)
# FIXME workaround for TP
monomial_type(v::Variables) = monomial_type(prod(v.variables))

constant_monomial_exponents(v::Variables) = map(_ -> 0, v.variables)

function (v::Variables)(exponents)
    return Polynomial(v, exponents)
end

variables(v::Variables) = v.variables
nvariables(v::Variables) = length(v.variables)

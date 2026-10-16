Variables{B}(vars) where {B} = Variables{B,typeof(vars)}(vars)

function Base.one(v::Variables)
    return Polynomial(v, constant_monomial_exponents(v))
end

function variable_index(v::Variables, var)
    return findfirst(isequal(var), v.variables)
end

function Base.:(==)(v::Variables{B}, w::Variables{B}) where {B}
    return v.variables === w.variables || v.variables == w.variables
end

function Base.hash(v::Variables{B}, u::UInt) where {B}
    return hash(v.variables, hash(B, u))
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

function _show(io::IO, mime::MIME, v::Variables{B}) where {B}
    print(io, "$B polynomials in the variables ")
    # We don't use the default `show` since we don't want to print the `eltype`
    # and we want to use the `mime`
    _show_vector(io, mime, v.variables)
    return
end

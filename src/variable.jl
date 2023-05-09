export name, name_base_indices, similar_variable, @similar_variable, variable_union_type

Base.copy(x::AbstractVariable) = x

"""
    variable_union_type(p::AbstractPolynomialLike)

Return the supertype for variables of `p`. If `p` is a variable, it should not
be the type of `p` but the supertype of all variables that could be created.

### Examples

For `TypedPolynomials`, a variable of name `x` has type `Variable{:x}` so
`variable_union_type` should return `Variable`.
For `DynamicPolynomials`, all variables have the same type `PolyVar{C}` where
`C` is `true` for commutative variables and `false` for non-commutative ones so
`variable_union_type` should return `PolyVar{C}`.
"""
function variable_union_type end
variable_union_type(::Type{MT}) where {MT<:AbstractMonomialLike} = error("`variable_union_type` not implemented for $MT.")
variable_union_type(::Type{PT}) where {PT<:APL} = variable_union_type(monomial_type(PT))
variable_union_type(p::APL) = variable_union_type(typeof(p))

"""
    variable(p::AbstractPolynomialLike)

Converts `p` to a variable. Throws `InexactError` if it is not possible.

### Examples

Calling `variable(x^2 + x - x^2)` should return the variable `x` and
calling `variable(1.0y)` should return the variable `y` however calling
`variable(2x)` or `variable(x + y)` should throw `InexactError`.

### Note

This operation is not type stable for the TypedPolynomials implementation if `nvariables(p) > 1` but is type stable for DynamicPolynomials.
"""
variable(t::APL) = convert(variable_union_type(t), t)

"""
    name(v::AbstractVariable)::AbstractString

Returns the name of a variable.
"""
function name end

"""
    name_base_indices(v::AbstractVariable)

Returns the name of the variable (as a `String` or `Symbol`) and its indices
as a `Vector{Int}` or tuple of `Int`s.
"""
function name_base_indices end

"""
    similar_variable(p::AbstractPolynomialLike, variable::Type{Val{V}})

Creates a new variable `V` based upon the the given source polynomial.

    similar_variable(p::AbstractPolynomialLike, v::Symbol)

Creates a new variable based upon the given source polynomial and the given symbol `v`. Note
that this can lead to type instabilities.

### Examples

Calling `similar_variable(typedpoly, Val{:x})` on a polynomial created with `TypedPolynomials`
results in `TypedPolynomials.Variable{:x}`.
"""
function similar_variable end

similar_variable(v::Union{AbstractVariable, Type{<:AbstractVariable}}, s::Symbol) = similar_variable(v, Val{s})
similar_variable(p::Union{AbstractPolynomialLike, Type{<:AbstractPolynomialLike}}, s) = similar_variable(variable_union_type(p), s)

"""
    @similar_variable(p::AbstractPolynomialLike, variable)

Calls `similar_variable(p, Val{variable})` and binds the result to a variable with the same
name.

### Examples

Calling `@similar_variable typedpoly x` on a polynomial created with `TypedPolynomials`
binds `TypedPolynomials.Variable{:x}` to the variable `x`.
"""
macro similar_variable(p, name)
    Expr(:block, _assign_new_variable(p, name))
end

_new_variable(p, name::Symbol) = :(similar_variable($p, Val{$(esc(Expr(:quote, name)))}))
_assign_new_variable(f, name::Symbol) = :($(esc(name)) = $(_new_variable(:($(esc(f))), name)))

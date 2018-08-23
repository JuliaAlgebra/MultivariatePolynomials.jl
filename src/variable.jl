export name, name_base_indices, similarvariable, @similarvariable

Base.copy(x::AbstractVariable) = x

"""
    variable(p::AbstractPolynomialLike)

Converts `p` to a variable. Throws an error if it is not possible.

### Examples

Calling `variable(x^2 + x - x^2)` should return the variable `x` and
calling `variable(1.0y)` should return the variable `y` however calling
`variable(2x)` or `variable(x + y)` should throw an error.

### Note

This operation is not type stable for the TypedPolynomials implementation if `nvariables(p) > 1` but is type stable for DynamicPolynomials.
"""
variable(m::AbstractMonomialLike) = _mono2var(powers(m)...)
variable(v::AbstractVariable) = v
function variable(t::AbstractTermLike)
    if isone(coefficient(t))
        variable(monomial(t))
    else
        error("A term with non-one coefficient cannot be converted into a variable")
    end
end
variable(p::APL) = variable(term(p))

_errormono2var() = error("Monomial cannot be converted to a variable")
_mono2var() = _errormono2var()
function _checknovar() end
function _checknovar(ve, ves...)
    if iszero(ve[2])
        _checknovar(ves...)
    else
        _errormono2var()
    end
end
function _mono2var(ve, ves...)
    if iszero(ve[2])
        _mono2var(ves...)
    elseif isone(ve[2])
        _checknovar(ves...)
        ve[1]
    else
        _errormono2var()
    end
end

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
    similarvariable(p::AbstractPolynomialLike, variable::Type{Val{V}})

Creates a new variable `V` based upon the the given source polynomial.

    similarvariable(p::AbstractPolynomialLike, v::Symbol)

Creates a new variable based upon the given source polynomial and the given symbol `v`. Note
that this can lead to type instabilities.

### Examples

Calling `similarvariable(typedpoly, Val{:x})` on a polynomial created with `TypedPolynomials`
results in `TypedPolynomials.Variable{:x}`.
"""
function similarvariable end

similarvariable(p::Union{AbstractPolynomialLike, Type{<:AbstractPolynomialLike}}, s::Symbol) = similarvariable(p, Val{s})

"""
    @similarvariable(p::AbstractPolynomialLike, variable)

Calls `similarvariable(p, Val{variable})` and binds the result to a variable with the same
name.

### Examples

Calling `@similarvariable typedpoly x` on a polynomial created with `TypedPolynomials`
binds `TypedPolynomials.Variable{:x}` to the variable `x`.
"""
macro similarvariable(p, name)
    Expr(:block, _makevar(p, name))
end

_varconstructor(p, name::Symbol) = :(similarvariable($p, Val{$(esc(Expr(:quote, name)))}))
_makevar(f, name::Symbol) = :($(esc(name)) = $(_varconstructor(:($(esc(f))), name)))

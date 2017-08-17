export similarvariable, @similarvariable

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

similarvariable(p::AbstractPolynomialLike, s::Symbol) = similarvariable(p, Val{s})


_varconstructor(p, name::Symbol) = :(similarvariable($p, Val{$(esc(Expr(:quote, name)))}))
_makevar(f, name::Symbol) = :($(esc(name)) = $(_varconstructor(:($(esc(f))), name)))

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

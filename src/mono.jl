export name, constantmonomial, monovec, monovectype, sortmonovec, mergemonovec

Base.copy(x::AbstractVariable) = x

"""
    name(v::AbstractVariable)::AbstractString

Returns the name of a variable.
"""
function name end

"""
    constantmonomial(p::AbstractPolynomialType)

Returns a constant monomial of the monomial type of `p` with the same variables as `p`.

    constantmonomial{PT<:AbstractPolynomialType}(::Type{PT})

Returns a constant monomial of the monomial type of a polynomial of type `PT`.
"""
function constantmonomial end

"""
    monovec{MT<:AbstractMonomialLike}(X::AbstractVector(MT)}

Returns the vector of monomials `X` in decreasing order and without any duplicates.

### Examples

Calling `monovec` on ``[xy, x, xy, x^2y, x]`` should return ``[x^2y, xy, x]``.
"""
function monovec end

"""
    monovectype{MT<:AbstractMonomialLike}(X::AbstractVector(MT)}

Returns the return type of `monovec`.
"""
function monovectype end

# If there are duplicates in X, the coefficients should be summed for a polynomial and they should be equal for a measure.
"""
    sortmonovec{MT<:AbstractMonomialLike}(X::AbstractVector(MT)}

Returns `σ`, the orders in which one must take the monomials in `X` to make them sorted and without any duplicate and the sorted vector of monomials, i.e. it returns `(σ, X[σ])`.

### Examples

Calling `sortmonovec` on ``[xy, x, xy, x^2y, x]`` should return ``([4, 1, 2], [x^2y, xy, x])``.
"""
function sortmonovec end

"""
    mergemonovec{MT<:AbstractMonomialLike, MVT<:AbstractVector{MT}}(X::Vector{MVT}}

Returns the vector of monomials in the entries of `X` in decreasing order and without any duplicates, i.e. `monovec(vat(X...))

### Examples

Calling `mergemonovec` on ``[[xy, x, xy], [x^2y, x]]`` should return ``[x^2y, xy, x]``.
"""
function mergemonovec end

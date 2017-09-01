export monovec, monovectype, emptymonovec, sortmonovec, mergemonovec


"""
    emptymonovec(p::AbstractPolynomialLike)

Returns an empty collection of the type of `monomials(p)`.
"""
emptymonovec(p) = monomialtype(p)[]

"""
    monovec(X::AbstractVector{MT}) where {MT<:AbstractMonomialLike}

Returns the vector of monomials `X` in decreasing order and without any duplicates.

### Examples

Calling `monovec` on ``[xy, x, xy, x^2y, x]`` should return ``[x^2y, xy, x]``.
"""
function monovec(X::AbstractVector{MT}) where {MT<:AbstractMonomial}
    Y = sort(X, rev=true)
    dups = find(i -> Y[i] == Y[i-1], 2:length(Y))
    deleteat!(Y, dups)
    Y
end
monovec(X::AbstractVector{TT}) where {TT<:AbstractTermLike} = monovec(AbstractVector{monomialtype(TT)}(X))

function monovec(a, x)
    if length(a) != length(x)
        throw(ArgumentError("There should be as many coefficient than monomials"))
    end
    σ, X = sortmonovec(x)
    (a[σ], X)
end

"""
    monovectype(X::AbstractVector{MT}) where {MT<:AbstractMonomialLike}

Returns the return type of `monovec`.
"""
monovectype(X::Union{AbstractVector{PT}, Type{<:AbstractVector{PT}}}) where {PT<:APL} = monovectype(PT)
monovectype(::Union{PT, Type{PT}}) where {PT <: APL} = Vector{monomialtype(PT)}

# If there are duplicates in X, the coefficients should be summed for a polynomial and they should be equal for a measure.
"""
    sortmonovec(X::AbstractVector{MT}) where {MT<:AbstractMonomialLike}

Returns `σ`, the orders in which one must take the monomials in `X` to make them sorted and without any duplicate and the sorted vector of monomials, i.e. it returns `(σ, X[σ])`.

### Examples

Calling `sortmonovec` on ``[xy, x, xy, x^2y, x]`` should return ``([4, 1, 2], [x^2y, xy, x])``.
"""
function sortmonovec(X::AbstractVector{MT}) where {MT<:AbstractMonomial}
    σ = sortperm(X, rev=true)
    dups = find(i -> X[σ[i]] == X[σ[i-1]], 2:length(σ))
    deleteat!(σ, dups)
    σ, X[σ]
end
sortmonovec(X::AbstractVector{TT}) where {TT<:AbstractTermLike} = sortmonovec(AbstractVector{monomialtype(TT)}(X))
sortmonovec(X::Tuple) = sortmonovec(vec(X))

"""
    mergemonovec{MT<:AbstractMonomialLike, MVT<:AbstractVector{MT}}(X::AbstractVector{MVT}}

Returns the vector of monomials in the entries of `X` in decreasing order and without any duplicates, i.e. `monovec(vcat(X...))`

### Examples

Calling `mergemonovec` on ``[[xy, x, xy], [x^2y, x]]`` should return ``[x^2y, xy, x]``.
"""
mergemonovec(X) = monovec(vcat(X...))

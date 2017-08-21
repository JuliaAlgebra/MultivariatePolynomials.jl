export name, constantmonomial, emptymonovec, monovec, monovectype, sortmonovec, mergemonovec, mapexponents

mapexponents(f, m1::AbstractMonomialLike, m2::AbstractMonomialLike) = mapexponents(f, monomial(m1), monomial(m2))

Base.copy(x::AbstractVariable) = x

"""
    name(v::AbstractVariable)::AbstractString

Returns the name of a variable.
"""
function name end

_hashpowers(u::UInt) = u
function _hashpowers(u::UInt, power::Tuple, powers...)
    if iszero(power[2])
        _hashpowers(u, powers...)
    else
        _hashpowers(hash(power, u), powers...)
    end
end
function Base.hash(m::AbstractMonomial, u::UInt)
    nnz = count(!iszero, exponents(m))
    if iszero(nnz)
        hash(1, u)
    elseif isone(nnz) && isone(degree(m))
        hash(variable(m), u)
    else
        _hashpowers(u, powers(m)...)
    end
end

Base.one(::Type{TT}) where {TT<:AbstractMonomialLike} = constantmonomial(TT)
Base.one(t::AbstractMonomialLike) = constantmonomial(t)

"""
    constantmonomial(p::AbstractPolynomialType)

Returns a constant monomial of the monomial type of `p` with the same variables as `p`.

    constantmonomial(::Type{PT}) where {PT<:AbstractPolynomialType}

Returns a constant monomial of the monomial type of a polynomial of type `PT`.
"""
function constantmonomial end

emptymonovec(::Type{PT}) where PT = monomialtype(PT)[]

"""
    monovec(X::AbstractVector{MT}) where {MT<:AbstractMonomialLike}

Returns the vector of monomials `X` in decreasing order and without any duplicates.

### Examples

Calling `monovec` on ``[xy, x, xy, x^2y, x]`` should return ``[x^2y, xy, x]``.
"""
function monovec(X::AbstractVector{MT}) where {MT<:AbstractMonomialLike}
    Y = sort(X, rev=true)
    dups = find(i -> Y[i] == Y[i-1], 2:length(Y))
    deleteat!(Y, dups)
    Y
end
monovec(X::AbstractVector{TT}) where {TT<:AbstractTerm} = monovec(AbstractVector{monomialtype(TT)}(X))

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
monovectype(X::AbstractVector{TT}) where {TT<:AbstractTermLike} = monovectype(TT)
monovectype(::Type{PT}) where {PT <: APL} = Vector{monomialtype(PT)}

# If there are duplicates in X, the coefficients should be summed for a polynomial and they should be equal for a measure.
"""
    sortmonovec(X::AbstractVector{MT}) where {MT<:AbstractMonomialLike}

Returns `σ`, the orders in which one must take the monomials in `X` to make them sorted and without any duplicate and the sorted vector of monomials, i.e. it returns `(σ, X[σ])`.

### Examples

Calling `sortmonovec` on ``[xy, x, xy, x^2y, x]`` should return ``([4, 1, 2], [x^2y, xy, x])``.
"""
function sortmonovec(X::AbstractVector{MT}) where {MT<:AbstractMonomialLike}
    σ = sortperm(X, rev=true)
    dups = find(i -> X[σ[i]] == X[σ[i-1]], 2:length(σ))
    deleteat!(σ, dups)
    σ, X[σ]
end
sortmonovec(X::AbstractVector{TT}) where {TT<:AbstractTerm} = sortmonovec(AbstractVector{monomialtype(TT)}(X))
sortmonovec(X::Tuple) = sortmonovec(vec(X))

"""
    mergemonovec{MT<:AbstractMonomialLike, MVT<:AbstractVector{MT}}(X::AbstractVector{MVT}}

Returns the vector of monomials in the entries of `X` in decreasing order and without any duplicates, i.e. `monovec(vcat(X...))`

### Examples

Calling `mergemonovec` on ``[[xy, x, xy], [x^2y, x]]`` should return ``[x^2y, xy, x]``.
"""
mergemonovec(X) = monovec(vcat(X...))

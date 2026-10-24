function monomials(v::AbstractVariable, degree, args...)
    return monomials(variables(v), degree, args...)
end

function monomials(vars::Tuple{Vararg{AbstractVariable}}, args...)
    return monomials(collect(vars), args...)
end

function monomials(
    vars::AbstractVector{<:AbstractVariable},
    degree::Int,
    args...,
)
    return monomials(vars, degree:degree, args...)
end

function monomials(
    vars::AbstractVector{<:AbstractVariable},
    degrees::AbstractVector{Int},
    filter::F = _ -> true,
) where {F}
    vars = unique!(sort(vars; rev = true))
    return _monomials(vars, degrees, filter)
end

function _monomials(vars, degrees, filter::F) where {F}
    if !is_commutative(vars)
        return _noncommutative_monomials(vars, degrees, filter)
    end
    basis = FullBasis{Monomial}(vars)
    mindeg, maxdeg = isempty(degrees) ? (0, -1) : extrema(degrees)
    exps = ExponentsIterator{typeof(ordering(vars))}(
        constant_monomial_exponents(basis);
        mindegree = mindeg,
        maxdegree = maxdeg,
    )
    monos = eltype(basis)[]
    for exp in exps
        sum(exp; init = 0) in degrees || continue
        mono = basis[exp]
        filter(mono) && push!(monos, mono)
    end
    return monos
end

function _noncommutative_monomials(vars, degrees, filter::F) where {F}
    if !(ordering(vars) isa Graded{LexOrder})
        throw(
            ArgumentError(
                "Noncommutative enumeration requires Graded{LexOrder}.",
            ),
        )
    end
    foreach(error_for_negative_degree, degrees)
    n = length(vars)
    maxdeg = maximum(degrees; init = 0)
    # Each new factor advances by at most n - 1 positions in the repeated
    # variable sequence; repeated adjacent factors stay at the same position.
    len = iszero(n) ? 0 : (n - 1) * maxdeg + 1
    repeated_vars = [vars[mod1(i, n)] for i in 1:len]
    exps = zeros(Int, len)
    monos = monomial_type(vars)[]
    for deg in unique!(sort!(collect(degrees)))
        _noncommutative_monomials!(
            monos,
            repeated_vars,
            exps,
            1,
            n,
            deg,
            filter,
        )
    end
    return monos
end

function _noncommutative_monomials!(
    monos,
    vars,
    exps,
    start,
    n,
    degree,
    filter::F,
) where {F}
    if iszero(degree)
        mono = monomial(vars, copy(exps))
        filter(mono) && push!(monos, mono)
    else
        # Visit later positions first to generate increasing lexicographic order.
        for i in (start+n-1):-1:start
            exps[i] += 1
            _noncommutative_monomials!(
                monos,
                vars,
                exps,
                i,
                n,
                degree - 1,
                filter,
            )
            exps[i] -= 1
        end
    end
    return
end

"""
    empty_monomial_vector(p::AbstractPolynomialLike)

Returns an empty collection of the type of `monomials(p)`.
"""
empty_monomial_vector(p) = monomial_type(p)[]

"""
    monomial_vector(X::AbstractVector{MT}) where {MT<:AbstractMonomialLike}

Returns the vector of monomials `X` in increasing order and without any duplicates.

### Examples

Calling `monomial_vector` on ``[xy, x, xy, x^2y, x]`` should return ``[x^2y, xy, x]``.
"""
function monomial_vector(X::AbstractVector{MT}) where {MT<:AbstractMonomial}
    Y = sort(X)
    dups = findall(i -> Y[i] == Y[i-1], 2:length(Y))
    deleteat!(Y, dups)
    return Y
end
function monomial_vector(X::AbstractVector{TT}) where {TT<:AbstractTermLike}
    return monomial_vector(AbstractVector{monomial_type(TT)}(X))
end

"""
    monomial_vector(a, X::AbstractVector{MT}) where {MT<:AbstractMonomialLike}

Returns `b, Y` where `Y` is the vector of monomials of `X` in increasing order
and without any duplicates and `b` is the vector of corresponding coefficients
in `a`, where coefficients of duplicate entries are summed together.

### Examples

Calling `monomial_vector` on ``[2, 1, 4, 3, -1], [xy, x, xy, x^2y, x]`` should return
``[3, 6, 0], [x^2y, xy, x]``.
"""
function monomial_vector(a, x)
    if length(a) != length(x)
        throw(
            ArgumentError("There should be as many coefficient than monomials"),
        )
    end
    σ, X = sort_monomial_vector(x)
    b = a[σ]
    if length(x) > length(X)
        rev = Dict(X[j] => j for j in eachindex(σ))
        for i in eachindex(x)
            j = rev[x[i]]
            if i != σ[j]
                b[j] += a[i]
            end
        end
    end
    return b, X
end

"""
    monomial_vector_type(X::AbstractVector{MT}) where {MT<:AbstractMonomialLike}

Returns the return type of `monomial_vector`.
"""
function monomial_vector_type(
    ::Union{AbstractVector{PT},Type{<:AbstractVector{PT}}},
) where {PT<:_APL}
    return monomial_vector_type(PT)
end
function monomial_vector_type(::Union{PT,Type{PT}}) where {PT<:_APL}
    return Vector{monomial_type(PT)}
end

# If there are duplicates in X, the coefficients should be summed for a polynomial and they should be equal for a measure.
"""
    sort_monomial_vector(X::AbstractVector{MT}) where {MT<:AbstractMonomialLike}

Returns `σ`, the orders in which one must take the monomials in `X` to make them sorted and without any duplicate and the sorted vector of monomials, i.e. it returns `(σ, X[σ])`.

### Examples

Calling `sort_monomial_vector` on ``[xy, x, xy, x^2y, x]`` should return ``([4, 1, 2], [x^2y, xy, x])``.
"""
function sort_monomial_vector(
    X::AbstractVector{MT},
) where {MT<:AbstractMonomial}
    σ = sortperm(X)
    dups = findall(i -> X[σ[i]] == X[σ[i-1]], 2:length(σ))
    deleteat!(σ, dups)
    return σ, X[σ]
end
function sort_monomial_vector(
    X::AbstractVector{TT},
) where {TT<:AbstractTermLike}
    return sort_monomial_vector(AbstractVector{monomial_type(TT)}(X))
end
sort_monomial_vector(X::Tuple) = sort_monomial_vector(vec(X))

"""
    merge_monomial_vectors{MT<:AbstractMonomialLike, MVT<:AbstractVector{MT}}(X::AbstractVector{MVT}}

Returns the vector of monomials in the entries of `X` in increasing order and without any duplicates, i.e. `monomial_vector(vcat(X...))`

### Examples

Calling `merge_monomial_vectors` on ``[[xy, x, xy], [x^2y, x]]`` should return ``[x^2y, xy, x]``.
"""
merge_monomial_vectors(X) = monomial_vector(reduce(vcat, X))

function error_for_negative_degree(deg)
    if deg < 0
        throw(
            ArgumentError(
                "The degree should be a nonnegative number but the provided degree `$deg` is negative.",
            ),
        )
    end
end

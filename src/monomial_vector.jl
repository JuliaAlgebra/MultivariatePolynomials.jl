export monomial_vector,
    monomial_vector_type,
    empty_monomial_vector,
    sort_monomial_vector,
    merge_monomial_vectors

export LexOrder, InverseLexOrder, Reverse, Graded

"""
    abstract type AbstractMonomialOrdering end

Abstract type for monomial ordering as defined in [CLO13, Definition 2.2.1, p. 55]

[CLO13] Cox, D., Little, J., & OShea, D.
*Ideals, varieties, and algorithms: an introduction to computational algebraic geometry and commutative algebra*.
Springer Science & Business Media, **2013**.
"""
abstract type AbstractMonomialOrdering end

"""
    compare(a, b, order::Type{<:AbstractMonomialOrdering})

Returns a negative number if `a < b`, a positive number if `a > b` and zero if `a == b`.
The comparison is done according to `order`.
"""
function compare end

"""
    struct LexOrder <: AbstractMonomialOrdering end

Lexicographic (Lex for short) Order often abbreviated as *lex* order as defined in [CLO13, Definition 2.2.3, p. 56]

The [`Graded`](@ref) version is often abbreviated as *grlex* order and is defined in [CLO13, Definition 2.2.5, p. 58]

[CLO13] Cox, D., Little, J., & OShea, D.
*Ideals, varieties, and algorithms: an introduction to computational algebraic geometry and commutative algebra*.
Springer Science & Business Media, **2013**.
"""
struct LexOrder <: AbstractMonomialOrdering end

"""
    struct InverseLexOrder <: AbstractMonomialOrdering end

Inverse Lex Order defined in [CLO13, Exercise 2.2.6, p. 61] where it is abbreviated as *invlex*.
It corresponds to [`LexOrder`](@ref) but with the variables in reverse order.

The [`Graded`](@ref) version can be abbreviated as *grinvlex* order.
It is defined in [BDD13, Definition 2.1] where it is called *Graded xel order*.

[CLO13] Cox, D., Little, J., & OShea, D.
*Ideals, varieties, and algorithms: an introduction to computational algebraic geometry and commutative algebra*.
Springer Science & Business Media, **2013**.
[BDD13] Batselier, K., Dreesen, P., & De Moor, B.
*The geometry of multivariate polynomial division and elimination*.
SIAM Journal on Matrix Analysis and Applications, 34(1), 102-125, *2013*.
"""
struct InverseLexOrder <: AbstractMonomialOrdering end

"""
    struct Graded{O<:AbstractMonomialOrdering} <: AbstractMonomialOrdering
        same_degree_ordering::O
    end

Monomial ordering defined by:
* `degree(a) == degree(b)` then the ordering is determined by `same_degree_ordering`,
* otherwise, it is the ordering between the integers `degree(a)` and `degree(b)`.
"""
struct Graded{O<:AbstractMonomialOrdering} <: AbstractMonomialOrdering
    same_degree_ordering::O
end

_deg(exponents) = sum(exponents)
_deg(mono::AbstractMonomial) = degree(mono)

function compare(a, b, ::Type{Graded{O}}) where {O}
    deg_a = _deg(a)
    deg_b = _deg(b)
    if deg_a == deg_b
        return compare(a, b, O)
    else
        return deg_a - deg_b
    end
end

"""
    struct Reverse{O<:AbstractMonomialOrdering} <: AbstractMonomialOrdering
        reverse_order::O
    end

Monomial ordering defined by
`compare(a, b, ::Type{Reverse{O}}) where {O} = compare(b, a, O)`.

Reverse Lex Order defined in [CLO13, Exercise 2.2.9, p. 61] where it is abbreviated as *rinvlex*.
can be obtained as `Reverse(InverseLexOrder())`.

The Graded Reverse Lex Order often abbreviated as *grevlex* order defined in [CLO13, Definition 2.2.6, p. 58]
can be obtained as `Graded(Reverse(InverseLexOrder()))`.

[CLO13] Cox, D., Little, J., & OShea, D.
*Ideals, varieties, and algorithms: an introduction to computational algebraic geometry and commutative algebra*.
Springer Science & Business Media, **2013**.
"""
struct Reverse{O<:AbstractMonomialOrdering} <: AbstractMonomialOrdering
    reverse_ordering::O
end

compare(a, b, ::Type{Reverse{O}}) where {O} = compare(b, a, O)

function monomials(v::AbstractVariable, degree, args...)
    return monomials((v,), degree, args...)
end

"""
    empty_monomial_vector(p::AbstractPolynomialLike)

Returns an empty collection of the type of `monomials(p)`.
"""
empty_monomial_vector(p) = monomial_type(p)[]

"""
    monomial_vector(X::AbstractVector{MT}) where {MT<:AbstractMonomialLike}

Returns the vector of monomials `X` in decreasing order and without any duplicates.

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

Returns `b, Y` where `Y` is the vector of monomials of `X` in decreasing order
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
    X::Union{AbstractVector{PT},Type{<:AbstractVector{PT}}},
) where {PT<:APL}
    return monomial_vector_type(PT)
end
function monomial_vector_type(::Union{PT,Type{PT}}) where {PT<:APL}
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

Returns the vector of monomials in the entries of `X` in decreasing order and without any duplicates, i.e. `monomial_vector(vcat(X...))`

### Examples

Calling `merge_monomial_vectors` on ``[[xy, x, xy], [x^2y, x]]`` should return ``[x^2y, xy, x]``.
"""
merge_monomial_vectors(X) = monomial_vector(reduce(vcat, X))

"""
    antidifferentiate(p::AbstractPolynomialLike, v::AbstractVariable, deg::Union{Int, Val}=1)

Antidifferentiate `deg` times the polynomial `p` by the variable `v`. The free constant involved by
the antidifferentiation is set to 0.


    antidifferentiate(p::AbstractPolynomialLike, vs, deg::Union{Int, Val}=1)

Antidifferentiate `deg` times the polynomial `p` by the variables of the vector or
tuple of variable `vs` and return an array of dimension `deg`. It is recommended
to pass `deg` as a `Val` instance when the degree is known at compile time, e.g.
`antidifferentiate(p, v, Val{2}())` instead of `antidifferentiate(p, x, 2)`, as this
will help the compiler infer the return type.

### Examples

```julia
p = 3x^2*y + x + 2y + 1
antidifferentiate(p, x) # should return 3x^3* + 1/2*x + 2xy + x
antidifferentiate(p, x, Val{1}()) # equivalent to the above
antidifferentiate(p, (x, y)) # should return [3x^3* + 1/2*x + 2xy + x, 3/2x^2*y^2 + xy + y^2 + y]
```
"""
function antidifferentiate end

# Fallback for everything else
antidifferentiate(α::T, v::AbstractVariable) where {T} = α * v
function antidifferentiate(v1::AbstractVariable, v2::AbstractVariable)
    return v1 == v2 ? 1 // 2 * v1 * v2 : v1 * v2
end
function antidifferentiate(t::AbstractTermLike, v::AbstractVariable)
    return coefficient(t) * antidifferentiate(monomial(t), v)
end
function antidifferentiate(p::_APL, v::AbstractVariable)
    return polynomial!(antidifferentiate.(terms(p), v), SortedState())
end

# TODO: this signature is probably too wide and creates the potential
# for stack overflows
antidifferentiate(p::_APL, xs) = [antidifferentiate(p, x) for x in xs]

# antidifferentiate(p, [x, y]) with TypedPolynomials promote x to a Monomial
function antidifferentiate(p::_APL, m::AbstractMonomial)
    return antidifferentiate(p, variable(m))
end

# The `R` argument indicates a desired result type. We use this in order
# to attempt to preserve type-stability even though the value of `deg` cannot
# be known at compile time. For scalar `p` and `x`, we set R to be the type
# of antidifferentiate(p, x) to give a stable result type regardless of `deg`. For
# vectors p and/or x this is impossible (since antidifferentiate may return an array),
# so we just set `R` to `Any`
function (_antidifferentiate_recursive(p, x, deg::Int, ::Type{R})::R) where {R}
    if deg < 0
        throw(DomainError(deg, "degree is negative"))
    elseif deg == 0
        return p
    else
        return antidifferentiate(antidifferentiate(p, x), x, deg - 1)
    end
end

function antidifferentiate(p, x, deg::Int)
    return _antidifferentiate_recursive(
        p,
        x,
        deg,
        Base.promote_op(antidifferentiate, typeof(p), typeof(x)),
    )
end
function antidifferentiate(p::AbstractArray, x, deg::Int)
    return _antidifferentiate_recursive(p, x, deg, Any)
end
function antidifferentiate(p, x::Union{AbstractArray,Tuple}, deg::Int)
    return _antidifferentiate_recursive(p, x, deg, Any)
end
function antidifferentiate(
    p::AbstractArray,
    x::Union{AbstractArray,Tuple},
    deg::Int,
)
    return _antidifferentiate_recursive(p, x, deg, Any)
end

# This is alternative, Val-based interface for nested antidifferentiation.
# It has the advantage of not requiring an conversion or calls to
# Base.promote_op, while maintaining type stability for any argument
# type.
antidifferentiate(p, x, ::Val{0}) = p
antidifferentiate(p, x, ::Val{1}) = antidifferentiate(p, x)

function antidifferentiate(p, x, deg::Val{N}) where {N}
    if N < 0
        throw(DomainError(deg))
    else
        antidifferentiate(antidifferentiate(p, x), x, Val{N - 1}())
    end
end

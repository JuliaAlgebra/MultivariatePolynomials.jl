export divides

"""
    divides(t1::AbstractTermLike, t2::AbstractTermLike)

Returns whether the monomial of t1 divides the monomial of t2.

### Examples

Calling `divides(2x^2y, 3xy)` should return false because `x^2y` does not divide `xy` since `x` has a degree 2 in `x^2y` which is greater than the degree of `x` on `xy`.
However, calling `divides(3xy, 2x^2y)` should return true.
"""
function divides(t1::AbstractTermLike, t2::AbstractTermLike)
    divides(monomial(t1), monomial(t2))
end
divides(t1::AbstractVariable, t2::AbstractVariable) = t1 == t2

Base.gcd(m1::AbstractMonomialLike, m2::AbstractMonomialLike) = mapexponents(min, m1, m2)
Base.lcm(m1::AbstractMonomialLike, m2::AbstractMonomialLike) = mapexponents(max, m1, m2)

# _div(a, b) assumes that b divides a
_div(m1::AbstractMonomialLike, m2::AbstractMonomialLike) = mapexponents(-, m1, m2)
function _div(t::AbstractTerm, m::AbstractMonomial)
    coefficient(t) * _div(monomial(t), m)
end
function _div(t1::AbstractTermLike, t2::AbstractTermLike)
    (coefficient(t1) / coefficient(t2)) * _div(monomial(t1), monomial(t2))
end

Base.div(f::APL, g::Union{APL, AbstractVector{<:APL}}; kwargs...) = divrem(f, g; kwargs...)[1]
Base.rem(f::APL, g::Union{APL, AbstractVector{<:APL}}; kwargs...) = divrem(f, g; kwargs...)[2]

function ring_rem(f::APL, g::APL)
    ltg = leadingterm(g)
    ltf = leadingterm(f)
    if !divides(monomial(ltg), ltf)
        return f
    end
    return coefficient(ltg) * removeleadingterm(f) -
           coefficient(ltf) * removeleadingterm(g)
end

proddiff(x, y) = x/y - x/y
function Base.divrem(f::APL{T}, g::APL{S}; kwargs...) where {T, S}
    # `promote_type(typeof(f), typeof(g))` is needed for TypedPolynomials in case they use different variables
    rf = convert(polynomialtype(promote_type(typeof(f), typeof(g)), Base.promote_op(proddiff, T, S)), MA.copy_if_mutable(f))
    q = zero(rf)
    r = zero(rf)
    lt = leadingterm(g)
    rg = removeleadingterm(g)
    lm = monomial(lt)
    while !iszero(rf)
        ltf = leadingterm(rf)
        if isapproxzero(ltf; kwargs...)
            rf = MA.operate!(removeleadingterm, rf)
        elseif divides(lm, ltf)
            qt = _div(ltf, lt)
            q = MA.add!(q, qt)
            rf = MA.operate!(removeleadingterm, rf)
            rf = MA.operate!(MA.sub_mul, rf, qt, rg)
        elseif lm > monomial(ltf)
            # Since the monomials are sorted in decreasing order,
            # lm is larger than all of them hence it cannot divide any of them
            r = MA.add!(r, rf)
            break
        else
            r = MA.add!(r, ltf)
            rf = MA.operate!(removeleadingterm, rf)
        end
    end
    q, r
end
function Base.divrem(f::APL{T}, g::AbstractVector{<:APL{S}}; kwargs...) where {T, S}
    # `promote_type(typeof(f), eltype(g))` is needed for TypedPolynomials in case they use different variables
    rf = convert(polynomialtype(promote_type(typeof(f), eltype(g)), Base.promote_op(proddiff, T, S)), MA.copy_if_mutable(f))
    r = zero(rf)
    q = similar(g, typeof(rf))
    for i in eachindex(q)
        q[i] = zero(rf)
    end
    lt = leadingterm.(g)
    rg = removeleadingterm.(g)
    lm = monomial.(lt)
    useful = BitSet(eachindex(g))
    while !iszero(rf)
        ltf = leadingterm(rf)
        if isapproxzero(ltf; kwargs...)
            rf = MA.operate!(removeleadingterm, rf)
            continue
        end
        divisionoccured = false
        for i in useful
            if divides(lm[i], ltf)
                qt = _div(ltf, lt[i])
                q[i] = MA.add!(q[i], qt)
                rf = MA.operate!(removeleadingterm, rf)
                rf = MA.operate!(MA.sub_mul, rf, qt, rg[i])
                divisionoccured = true
                break
            elseif lm[i] > monomial(ltf)
                # Since the monomials are sorted in decreasing order,
                # lm is larger than all of them hence it cannot divide any of them
                delete!(useful, i)
            end
        end
        if !divisionoccured
            if isempty(useful)
                r = MA.add!(r, rf)
                break
            else
                r = MA.add!(r, ltf)
                rf = MA.operate!(removeleadingterm, rf)
            end
        end
    end
    q, r
end

function isolate_variable(poly::APL, var::AbstractVariable)
    T = typeof(substitute(Subs(), zeroterm(poly), (var,) => (1,)))
    dict = Dict{Int,Vector{T}}()
    for t in terms(poly)
        d = degree(t, var)
        dict_ts = get(dict, d, nothing)
        ts = dict_ts === nothing ? T[] : dict_ts
        push!(ts, subs(t, (var,) => (1,)))
        if dict_ts === nothing
            dict[d] = ts
        end
    end
    return polynomial([
        term(polynomial(ts), var^d) for (d, ts) in dict
    ])
end

function coefficients_gcd(poly::APL{T}) where {T}
    return reduce(gcd, coefficients(poly), init=zero(T))
end

"""
    function Base.gcd(p1::AbstractPolynomialLike{T}, p2::AbstractPolynomialLike{S}) where {T, S}

Returns the greatest common divisor of `p1` and `p2`.

# Implementation notes

The classical algorithm for computing the `gcd`, commonly referred to as the
Euclidean Algorithm is to use a recursion with the base case `gcd(p, 0) = p`
and the relation `gcd(p1, p2) = gcd(p2, rem(p1, p2))`.
The relation comes from the Euclidean division:
`p1 = q * p2 + r`,
if `g` divides `p1` and `p2` then it divides `r` and
if `g` divides `r` and `p2` then it divides `p1`.

For multivariate polynomials, you may have `rem(p1, p2) = p1` hence this will not
terminate.
To ensure we make progress, we can pick a given variable `xi` and try to find
`q1` and `q2` such that
`q2 * p1 = q1 * p2 + r`
and the degree of `r` in `xi` is strictly smaller than the degree of `p1` in `xi`.
Note that if `g` divides `p1` and `p2` then it divides `r` but
if `g` divides `r` and `p2` then it might divide `q2` and not `p1`.
So what do we do ?
Let `dj` be the degree of `pj` in `xi`.
Suppose we pick `qj` to be the coefficient of `pj` in `xi^dj`.
If `g` divides `q2` then it means that the degree of `g` in `xi` is zero.
Therefore, if it divides `p2` then it also divides the coefficients
of `p2` in `xi^k` for `k = 0, 1, ..., d2`.
This means that if we ensure that these are relatively prime then we won't have
any issue.
So we start by computing the `gcd` `gj` of the coefficients in each degree of
`xi` of `pj`.
And then we compute `_gcd(p1 / g1, p2 / g2) * gcd(g1, g2)` where we can use the
recursion `_gcd(p1, p2) = _gcd(p2, q2 * p1 - q1 * p2)` where `q1, q2` are as
defined above. The base case is that `_gcd(p1, p2)` redirects to `gcd(p1, p2)`
when the degree of both `p1` and `p2` is zero for `xi`. Then we pick another
variable and we continue until we arrive at `gcd(p, 0)`.
"""
function Base.gcd(p::AbstractPolynomialLike{T}, q::AbstractPolynomialLike{S}) where {T, S}
    p_vars = effective_variables(p)
    q_vars = effective_variables(q)
    if length(p_vars) == length(q_vars) == 1 && first(p_vars) == first(q_vars)
        return univariate_gcd(p, q)
    else
        vars = _first_in_union_decreasing(p_vars, q_vars)
        if isempty(vars)
            return univariate_gcd(p, q)
            #return one(MA.promote_operation(gcd, typeof(p), typeof(q)))
        else
            return multivariate_gcd(p, q, first(vars))
        end
    end
end
# Returns first element in the union of two decreasing vectors
function _first_in_union_decreasing(a, b)
    i = j = 1
    while i <= length(a) && j <= length(b)
        if a[i] == b[j]
            return a[i]
        elseif a[i] > b[j]
            i += 1
        else
            j += 1
        end
    end
    return
end

function MA.promote_operation(::typeof(gcd), P::Type{<:AbstractPolynomialLike{T}}, ::Type{<:AbstractPolynomialLike{S}}) where {T,S}
    return polynomialtype(P, Base.promote_op(proddiff, T, S))
end

function univariate_gcd(p::AbstractPolynomialLike, q::AbstractPolynomialLike)
    if isapproxzero(q)
        convert(MA.promote_operation(gcd, typeof(p), typeof(q)), p)
    else
        univariate_gcd(q, rem(p, q))
    end
end
function extracted_multivariate_gcd(p::AbstractPolynomialLike, q::AbstractPolynomialLike)
    if isapproxzero(q)
        return p
    else
        return extracted_multivariate_gcd(q, ring_rem(p, q))
    end
end
function multivariate_gcd(p1::AbstractPolynomialLike, p2::AbstractPolynomialLike, var)
    e1 = isolate_variable(p1, var)
    g1 = coefficients_gcd(e1)
    f1 = mapcoefficientsnz(e1) do α
        div(α, g1)
    end
    e2 = isolate_variable(p2, var)
    g2 = coefficients_gcd(e2)
    f2 = mapcoefficientsnz(e2) do α
        div(α, g2)
    end
    ee = extracted_multivariate_gcd(f1, f2)
    return sum(coefficient(t) * monomial(t) for t in terms(ee)) * gcd(g1, g2)
end

Base.lcm(p::AbstractPolynomialLike, q::AbstractPolynomialLike) = p * div(q, gcd(p, q))

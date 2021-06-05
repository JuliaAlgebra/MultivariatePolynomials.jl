Base.lcm(p::APL, q::APL) = p * div(q, gcd(p, q))
Base.gcd(α, p::APL) = gcd(α, content(p))
Base.gcd(p::APL, α) = gcd(content(p), α)

#function MA.promote_operation(::typeof(gcd), P::Type{<:Number}, Q::Type{<:Number})
#    return typeof(gcd(one(P), one(Q)))
#end
function MA.promote_operation(::typeof(gcd), P::Type{<:APL}, Q::Type{<:APL})
    return MA.promote_operation(ring_rem, P, Q)
end

"""
    function gcd(p1::AbstractPolynomialLike{T}, p2::AbstractPolynomialLike{S}) where {T, S}

Returns a greatest common divisor of `p1` and `p2`. Note that it does not make
sense, in general, to speak of "the" greatest common divisor of u and v; there
is a set of greatest common divisors, each one being a unit multiple of the
others [Knu14, p. 424].

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
So we start by computing a `gcd` `gj` of the coefficients in each degree of
`xi` of `pj`, this is called the [`content`](@ref) of `pj`.
And then we compute `_gcd(p1 / g1, p2 / g2) * gcd(g1, g2)` where we can use the
recursion `_gcd(p1, p2) = _gcd(p2, q2 * p1 - q1 * p2)` where `q1, q2` are as
defined above.
This is the [`GeneralizedEuclideanAlgorithm`](@ref).

[Knu14] Knuth, D.E., 2014.
*Art of computer programming, volume 2: Seminumerical algorithms.*
Addison-Wesley Professional. Third edition.
"""
function Base.gcd(p1::APL{T}, p2::APL{S}) where {T, S}
    # If one of these is zero, `shift` should be infinite
    # for this method to work so we exclude these cases.
    if isapproxzero(p1)
        return convert(MA.promote_operation(gcd, typeof(p1), typeof(p2)), p2)
    end
    if isapproxzero(p1)
        return convert(MA.promote_operation(gcd, typeof(p1), typeof(p2)), p2)
    end
    shift1, defl1 = deflation(p1)
    shift2, defl2 = deflation(p2)
    shift = gcd(shift1, shift2)
    defl = mapexponents(gcd, defl1, defl2)
    # We factor out `x.^shift1` from `p1` and
    # `x.^shift2` from `p2`. The `gcd` of these
    # monomials is `x.^shift`.
    # Then, we subsitute `y[i] = x[i]^defl[i]`.
    q1 = deflate(p1, shift1, defl)
    q2 = deflate(p2, shift2, defl)
    g = deflated_gcd(q1, q2)
    return inflate(g, shift, defl)::MA.promote_operation(gcd, typeof(p1), typeof(p2))
end

# Inspired from to `AbstractAlgebra.deflation`
function deflation(p::AbstractPolynomialLike)
    if iszero(p)
        return constantmonomial(p), constantmonomial(p)
    end
    shift = fill(-1, nvariables(p))
    defl = zeros(Int, nvariables(p))
    for mono in monomials(p)
        exps = exponents(mono)
        for i in eachindex(exps)
            exp = exps[i]
            @assert exp >= 0
            if shift[i] == -1
                shift[i] = exp
            elseif exp < shift[i]
                # There are two cases:
                # 1) If `defl[i]` is zero then it means all previous monomials
                #    had degree `shift[i]` so we just set `defl[i]` to
                #    `shift[i] - exp` or equivalently `gcd(0, shift[i] - exp)`.
                # 2) If  `defl[i]` is positive then we have some monomials with
                #    degree `shift[i]` and some with degree
                #    `shift[i] + k * defl[i]` for some `k > 0`. We have
                #    `gcd(shift[i] - exp, shift[i] + k1 * defl[i] - exp, shift[i] + k2 * defl[i] - exp, ...) =`
                #    `gcd(shift[i] - exp, k1 * defl[i], k2 * defl[i], ...)`
                #    Since `gcd(k1, k2, ...) = 1`, this is equal to
                #    `gcd(shift[i] - exp, defl[i])`
                defl[i] = gcd(defl[i], shift[i] - exp)
                shift[i] = exp
            else
                defl[i] = gcd(defl[i], exp - shift[i])
            end
        end
    end
    @assert all(d -> d >= 0, shift)
    @assert all(d -> d >= 0, defl)
    s = prod(variables(p).^shift)::monomialtype(p)
    d = prod(variables(p).^defl)::monomialtype(p)
    return s, d
end

function deflate(p::AbstractPolynomialLike, shift, defl)
    if isconstant(shift) && all(d -> isone(d) || iszero(d), exponents(defl))
        return p
    end
    q = MA.operate(deflate, p, shift, defl)
    return q
end
function inflate(α, shift, defl)
    return inflate(convert(polynomialtype(shift, typeof(α)), α), shift, defl)
end
function inflate(p::AbstractPolynomialLike, shift, defl)
    if isconstant(shift) && all(d -> isone(d) || iszero(d), exponents(defl))
        return p
    end
    q = MA.operate(inflate, p, shift, defl)
    return q
end

function MA.operate(::typeof(deflate), mono::AbstractMonomial, shift, defl)
    mutable_mono = mapexponents(-, mono, shift)
    return mapexponents!(div, mutable_mono, defl)
end
function MA.operate(::typeof(inflate), mono::AbstractMonomial, shift, defl)
    mutable_mono = mapexponents(*, mono, defl)
    return mapexponents!(+, mutable_mono, shift)
end

# Inspired from to `AbstractAlgebra.deflate`
function MA.operate(op::Union{typeof(deflate), typeof(inflate)}, p::AbstractPolynomialLike, shift, defl)
    defl = prod(variables(defl).^map(d -> iszero(d) ? one(d) : d, exponents(defl)))
    return polynomial(map(terms(p)) do t
        return term(coefficient(t), MA.operate(op, monomial(t), shift, defl))
    end)
end

function deflated_gcd(p1::APL{T}, p2::APL{S}) where {T, S}
    i1, i2, num_common = _extracted_variable(p1, p2)
    if iszero(i1)
        if iszero(i2)
            return univariate_gcd(p1, p2)
        else
            if isapproxzero(p1)
                return convert(MA.promote_operation(gcd, typeof(p1), typeof(p2)), p2)
            end
            v2 = variables(p2)[i2]
            q2 = isolate_variable(p2, v2)
            g = content(q2)
            return gcd(p1, g)
        end
    else
        if iszero(i2)
            if isapproxzero(p2)
                return convert(MA.promote_operation(gcd, typeof(p1), typeof(p2)), p1)
            end
            v1 = variables(p1)[i1]
            q1 = isolate_variable(p1, v1)
            g = content(q1)
            return gcd(g, p2)
        else
            if num_common > 1
                @assert i1 == i2
                v1 = variables(p1)[i1]
                return multivariate_gcd(p1, p2, v1)
            else
                return univariate_gcd(p1, p2)
            end
        end
    end
end

# Returns first element in the union of two decreasing vectors
function _extracted_variable(p1, p2)
    v1 = variables(p1)
    v2 = variables(p2)
    i1 = i2 = 1
    best = nothing
    best_var1 = 0
    best_var2 = 0
    num_common = 0
    while i1 <= length(v1) || i2 <= length(v2)
        if i2 > length(v2) || (i1 <= length(v1) && v1[i1] > v2[i2])
            if !iszero(maxdegree(p1, v1[i1]))
                return i1, 0, num_common
            end
            i1 += 1
        elseif i1 > length(v1) || v2[i2] > v1[i1]
            if !iszero(maxdegree(p2, v2[i2]))
                return 0, i2, num_common
            end
            i2 += 1
        else
            @assert v1[i1] == v2[i2]
            v = v1[i1]
            i1 += 1
            i2 += 1
            d1, n1 = deg_num_leading_terms(p1, v)
            d2, n2 = deg_num_leading_terms(p2, v)
            if iszero(d1)
                if iszero(d2)
                    continue
                else
                    return 0, i2-1, num_common
                end
            else
                if iszero(d2)
                    return i1-1, 0, num_common
                end
            end
            if d1 < d2
                d1, d2 = d2, d1
                n1, n2 = n2, n1
            end
            # Heuristic used in `AbstractAlgebra`:
            # https://github.com/Nemocas/AbstractAlgebra.jl/blob/4c6b0a366e550df3db84a665de186111bc3cf8ed/src/generic/MPoly.jl#L4347
            # FIXME what is this based on ? Is there any analysis somewhere comparing different heuristics ?
            cur = max(log(n2) * d1 * d2, log(2) * d2)
            if best === nothing || best > cur
                best = cur
                best_var1 = i1-1
                best_var2 = i2-1
            end
            num_common += 1
        end
    end
    return best_var1, best_var2, num_common
end

function multivariate_gcd(p1::APL, p2::APL, var)
    q = univariate_gcd(isolate_variable(p1, var), isolate_variable(p2, var))
    return sum(coefficient(t) * monomial(t) for t in terms(q))::MA.promote_operation(gcd, typeof(p1), typeof(p2))
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

"""
    abstract type AbstractUnivariateGCDAlgorithm end

Algorithm computing the greatest common divisor of univariate polynomials.
See [`GeneralizedEuclideanAlgorithm`](@ref) and [`SubresultantAlgorithm`](@ref).
"""
abstract type AbstractUnivariateGCDAlgorithm end

"""
    struct GeneralizedEuclideanAlgorithm <: AbstractUnivariateGCDAlgorithm end

Algorithm computing the greatest common divisor of univariate polynomials using
the Euclidean algorithm generalized for polynomials with coefficients over a
a unique factorization domain, see [Knu14, Algorithm E, p. 426-427].

[Knu14] Knuth, D.E., 2014.
*Art of computer programming, volume 2: Seminumerical algorithms.*
Addison-Wesley Professional. Third edition.
"""
struct GeneralizedEuclideanAlgorithm <: AbstractUnivariateGCDAlgorithm end

"""
    struct SubresultantAlgorithm <: AbstractUnivariateGCDAlgorithm end

Algorithm computing the greatest common divisor of univariate polynomials using
the Subresultant algorithm, see [Knu14, Algorithm C, p. 428-429].

[Knu14] Knuth, D.E., 2014.
*Art of computer programming, volume 2: Seminumerical algorithms.*
Addison-Wesley Professional. Third edition.
"""
struct SubresultantAlgorithm <: AbstractUnivariateGCDAlgorithm end

"""
    primitive_univariate_gcd(p::APL, q::APL, algo::AbstractUnivariateGCDAlgorithm)

Returns the `gcd` of primitive polynomials `p` and `q` using algorithm `algo`
which is a subtype of [`AbstractUnivariateGCDAlgorithm`](@ref).
"""
function primitive_univariate_gcd end

function primitive_univariate_gcd(p::APL, q::APL, ::GeneralizedEuclideanAlgorithm=GeneralizedEuclideanAlgorithm())
    R = MA.promote_operation(gcd, typeof(p), typeof(q))
    if isapproxzero(q)
        return convert(R, p)
    elseif isapproxzero(p)
        return convert(R, q)
    elseif isconstant(p) || isconstant(q)
        # `p` and `q` are primitive so if one of them is constant, it cannot
        # divide the content of the other one.
        return one(R)
    else
        divided, r = ring_rem(p, q)
        o = q
        if !divided
            divided, r = ring_rem(q, p)
            o = p
            # Since `p` and `q` are univariate, at least one divides the other
            @assert divided
        end
        return primitive_univariate_gcd(o, primitive_part(r))
    end
end

function primitive_univariate_gcd(p::APL, q::APL, ::SubresultantAlgorithm)
    error("Not implemented yet")
end

"""
    univariate_gcd(p1::AbstractPolynomialLike, p2::AbstractPolynomialLike)

Return the *greatest common divisor* of the polynomials `p1` and `p2` that have
at most one variable in common and for which the coefficients are either
`AbstractFloat` or part of a unique factorization domain, e.g., rational numbers,
integers or multivariate polynomials. So `p1` and `p2` should have at most one
variable in common but their coefficients can be multivariate polynomials that
share arbitrarily many variables.

If the coefficients are not `AbstractFloat`, this
1. separates `p1` and `p2` in their [`content`](@ref) and
   [`primitive_part`](@ref) using [`primitive_part_content`](@ref); see
   [Knu14, Algorithm E: E1, p. 426] or [Knu14, Algorithm C: C1, p. 428].
2. Computes the [`gcd`](@ref) of the contents and primitive parts, using
   [`primitive_univariate_gcd`](@ref) for primitive parts.
3. Return the product of these two `gcd`; see
   [Knu14, Algorithm E: E4, p. 427] or [Knu14, Algorithm C: C4, p. 429].

[Knu14] Knuth, D.E., 2014.
*Art of computer programming, volume 2: Seminumerical algorithms.*
Addison-Wesley Professional. Third edition.
"""
function univariate_gcd(p1::APL, p2::APL)
    f1, g1 = primitive_part_content(p1)
    f2, g2 = primitive_part_content(p2)
    pp = primitive_univariate_gcd(f1, f2)
    gg = gcd(g1, g2)#::MA.promote_operation(gcd, typeof(g1), typeof(g2))
    # Multiply each coefficient by the gcd of the contents.
    return mapcoefficientsnz(Base.Fix1(*, gg), pp)
end

function univariate_gcd(p1::APL, p2::APL{<:AbstractFloat})
    return primitive_univariate_gcd(p1, p2)
end

_simplifier(a, b) = gcd(a, b)
# Before Julia v1.4, it is not defined.
# After Julia v1.4, it is defined as `gcd of num / lcm of den`.
# We prefer `gcd of den`, otherwise,
# `1/a0 + 1/a1 x + 1/a2 x^2 + ... + 1/an x^n`
# will be transformed into
# `a1*a2*...*an + a0*a2*...*an x + ...`
# which makes the size of the `BigInt`s grow significantly which slows things down.
_simplifier(a::Rational, b::Rational) = gcd(a.num, b.num) // gcd(a.den, b.den)

"""
    content(poly::AbstractPolynomialLike{T}) where {T}

Return the *content* of the polynomial `poly` over a unique factorization
domain `S` as defined in [Knu14, (3) p. 423].
That is, return the `gcd` of the coefficients of `poly`.
See also [`primitive_part_content`](@ref).

[Knu14] Knuth, D.E., 2014.
*Art of computer programming, volume 2: Seminumerical algorithms.*
Addison-Wesley Professional. Third edition.
"""
function content(poly::APL{T}) where {T}
    P = MA.promote_operation(gcd, T, T)
    # This is tricky to infer a `content` calls `gcd` which calls `content`, etc...
    # To help Julia break the loop, we annotate the result here.
    return reduce(
        _simplifier,
        coefficients(poly),
        init = zero(P),
    )::P
end

"""
    primitive_part_content(poly::AbstractPolynomialLike{T}) where {T}

Return the *primitive part* of the polynomial `poly` over a unique
factorization domain `S` as defined in [Knu14, (3) p. 423].
That is, return the exact division of `poly` by its [`content`](@ref).
If the content is also needed, call [`primitive_part_content`](@ref)
instead.

[Knu14] Knuth, D.E., 2014.
*Art of computer programming, volume 2: Seminumerical algorithms.*
Addison-Wesley Professional. Third edition.
"""
primitive_part(p::APL) = primitive_part_content(p)[1]
primitive_part(p::APL{<:AbstractFloat}) = p

"""
    primitive_part_content(poly::AbstractPolynomialLike{T}) where {T}

Return the *primitive part* and *content* of the polynomial `poly` over a unique
factorization domain `S` as defined in [Knu14, (3) p. 423]. This is more
efficient to call this function rather than calling [`primitive_part`](@ref) and
[`content`](@ref) separately since computing the primitive part requires
computing the content first and this function avoid computing the content twice.

[Knu14] Knuth, D.E., 2014.
*Art of computer programming, volume 2: Seminumerical algorithms.*
Addison-Wesley Professional. Third edition.
"""
function primitive_part_content(p)
    g = content(p)
    return mapcoefficientsnz(Base.Fix2(_div, g), p), g
end

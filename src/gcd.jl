export GeneralizedEuclideanAlgorithm

"""
    abstract type AbstractUnivariateGCDAlgorithm end

Algorithm computing the greatest common divisor of univariate polynomials.
See [`GeneralizedEuclideanAlgorithm`](@ref) and [`SubresultantAlgorithm`](@ref).
"""
abstract type AbstractUnivariateGCDAlgorithm end

"""
    struct GeneralizedEuclideanAlgorithm <: AbstractUnivariateGCDAlgorithm
        primitive_rem::Bool
        skip_last::Bool
    end

Algorithm computing the greatest common divisor of univariate polynomials using
the Euclidean algorithm generalized for polynomials with coefficients over a
a unique factorization domain, see [Knu14, Algorithm E, p. 426-427].

[Knu14] Knuth, D.E., 2014.
*Art of computer programming, volume 2: Seminumerical algorithms.*
Addison-Wesley Professional. Third edition.
"""
struct GeneralizedEuclideanAlgorithm <: AbstractUnivariateGCDAlgorithm
    primitive_rem::Bool
    skip_last::Bool
    function GeneralizedEuclideanAlgorithm(
        primitive_rem::Bool=false,
        skip_last::Bool=false,
    )
        return new(primitive_rem, skip_last)
    end
end

"""
    struct SubresultantAlgorithm <: AbstractUnivariateGCDAlgorithm end

Algorithm computing the greatest common divisor of univariate polynomials using
the Subresultant algorithm, see [Knu14, Algorithm C, p. 428-429].

[Knu14] Knuth, D.E., 2014.
*Art of computer programming, volume 2: Seminumerical algorithms.*
Addison-Wesley Professional. Third edition.
"""
struct SubresultantAlgorithm <: AbstractUnivariateGCDAlgorithm end

_coefficient_gcd(α, β) = gcd(α, β)
_coefficient_gcd(α::AbstractFloat, β) = one(Base.promote_typeof(α, β))
_coefficient_gcd(α, β::AbstractFloat) = one(Base.promote_typeof(α, β))
_coefficient_gcd(α::AbstractFloat, β::AbstractFloat) = one(Base.promote_typeof(α, β))

Base.lcm(p::APL, q::APL, algo::AbstractUnivariateGCDAlgorithm=GeneralizedEuclideanAlgorithm()) = p * div(q, gcd(p, q, algo))
Base.gcd(α, p::APL, algo::AbstractUnivariateGCDAlgorithm=GeneralizedEuclideanAlgorithm()) = _coefficient_gcd(α, content(p, algo))
Base.gcd(p::APL, α, algo::AbstractUnivariateGCDAlgorithm=GeneralizedEuclideanAlgorithm()) = _coefficient_gcd(content(p, algo), α)

#function MA.promote_operation(::typeof(gcd), P::Type{<:Number}, Q::Type{<:Number})
#    return typeof(gcd(one(P), one(Q)))
#end
function MA.promote_operation(::typeof(gcd), P::Type{<:APL}, Q::Type{<:APL})
    return MA.promote_operation(pseudo_rem, P, Q)
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
function Base.gcd(p1::APL{T}, p2::APL{S}, algo::AbstractUnivariateGCDAlgorithm=GeneralizedEuclideanAlgorithm()) where {T, S}
    # If one of these is zero, `shift` should be infinite
    # for this method to work so we exclude these cases.
    if isapproxzero(p1)
        return convert(MA.promote_operation(gcd, typeof(p1), typeof(p2)), p2)
    end
    if isapproxzero(p2)
        return convert(MA.promote_operation(gcd, typeof(p1), typeof(p2)), p1)
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
    g = deflated_gcd(q1, q2, algo)
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

function _zero_to_one_exp(defl::AbstractMonomial)
    # TODO Make it faster by calling something like `mapexponents`.
    return prod(variables(defl).^map(d -> iszero(d) ? one(d) : d, exponents(defl)))
end
function deflate(p::AbstractPolynomialLike, shift, defl)
    if isconstant(shift) && all(d -> isone(d) || iszero(d), exponents(defl))
        return p
    end
    q = MA.operate(deflate, p, shift, _zero_to_one_exp(defl))
    return q
end
function inflate(α, shift, defl)
    return inflate(convert(polynomialtype(shift, typeof(α)), α), shift, defl)
end
function inflate(p::AbstractPolynomialLike, shift, defl)
    if isconstant(shift) && all(d -> isone(d) || iszero(d), exponents(defl))
        return p
    end
    q = MA.operate(inflate, p, shift, _zero_to_one_exp(defl))
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
    return polynomial(map(terms(p)) do t
        return term(coefficient(t), MA.operate(op, monomial(t), shift, defl))
    end)
end

function deflated_gcd(p1::APL{T}, p2::APL{S}, algo) where {T, S}
    i1, i2, num_common = _extracted_variable(p1, p2)
    if iszero(i1)
        if iszero(i2)
            return univariate_gcd(p1, p2, algo)
        else
            if isapproxzero(p1)
                return convert(MA.promote_operation(gcd, typeof(p1), typeof(p2)), p2)
            end
            v2 = variables(p2)[i2]
            q2 = isolate_variable(p2, v2)
            g = content(q2, algo)
            return gcd(p1, g, algo)
        end
    else
        if iszero(i2)
            if isapproxzero(p2)
                return convert(MA.promote_operation(gcd, typeof(p1), typeof(p2)), p1)
            end
            v1 = variables(p1)[i1]
            q1 = isolate_variable(p1, v1)
            g = content(q1, algo)
            return gcd(g, p2, algo)
        else
            if num_common > 1
                v1 = variables(p1)[i1]
                @assert v1 == variables(p2)[i2]
                return multivariate_gcd(p1, p2, v1, algo)
            else
                return univariate_gcd(p1, p2, algo)
            end
        end
    end
end

function Base.gcdx(p1::APL{T}, p2::APL{S}, algo::AbstractUnivariateGCDAlgorithm=GeneralizedEuclideanAlgorithm()) where {T, S}
    i1, i2, num_common = _extracted_variable(p1, p2)
    R = MA.promote_operation(gcd, typeof(p1), typeof(p2))
    if iszero(i1)
        if iszero(i2)
            return univariate_gcdx(p1, p2, algo)
        else
            if isapproxzero(p1)
                return zero(R), one(R), convert(R, p2)
            end
            error("Not implemented yet")
        end
    else
        if iszero(i2)
            if isapproxzero(p2)
                return one(R), zero(R), convert(R, p1)
            end
            error("Not implemented yet")
        else
            if num_common > 1
                @assert i1 == i2
                error("Not implemented yet")
            else
                return univariate_gcdx(p1, p2, algo)
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

function multivariate_gcd(p1::APL, p2::APL, var, algo)
    q = univariate_gcd(isolate_variable(p1, var), isolate_variable(p2, var), algo)
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
    primitive_univariate_gcd(p::APL, q::APL, algo::AbstractUnivariateGCDAlgorithm)

Returns the `gcd` of primitive polynomials `p` and `q` using algorithm `algo`
which is a subtype of [`AbstractUnivariateGCDAlgorithm`](@ref).
"""
function primitive_univariate_gcd end

function not_divided_error(u, v)
    error(
        "Polynomial `$v` of degree `$(maxdegree(v))` and effective",
        " variables `$(effective_variables(v))` does not divide",
        " polynomial `$u` of degree `$(maxdegree(u))` and effective",
        " variables `$(effective_variables(u))`. Did you call",
        " `univariate_gcd` with polynomials with more than one",
        " variable in common ? If yes, call `gcd` instead, otherwise,",
        " please report this.",
    )
end

# If `p` and `q` do not have the same time then the local variables `p` and `q`
# won't be type stable.
function primitive_univariate_gcd(p::APL, q::APL, algo::GeneralizedEuclideanAlgorithm)
    if maxdegree(p) < maxdegree(q)
        return primitive_univariate_gcd(q, p, algo)
    end
    R = MA.promote_operation(gcd, typeof(p), typeof(q))
    u = convert(R, p)
    v = convert(R, q)
    while true
        if isapproxzero(v)
            return u
        elseif isconstant(v)
            # `p` and `q` are primitive so if one of them is constant, it cannot
            # divide the content of the other one.
            return one(R)
        end
        divided, r = pseudo_rem(u, v, algo)
        if !divided
            not_divided_error(u, v)
        end
        if !algo.primitive_rem
            r = primitive_part(r, algo)::R
        end
        u, v = v, r::R
    end
end

function primitive_univariate_gcdx(u0::APL, v0::APL, algo::GeneralizedEuclideanAlgorithm)
    if maxdegree(u0) < maxdegree(v0)
        a, b, g = primitive_univariate_gcdx(v0, u0, algo)
        return b, a, g
    end
    R = MA.promote_operation(gcd, typeof(u0), typeof(v0))
    u = convert(R, u0)
    v = convert(R, v0)
    if isapproxzero(v)
        return one(R), zero(R), u
    elseif isconstant(v)
        # `p` and `q` are primitive so if one of them is constant, it cannot
        # divide the content of the other one.
        return zero(R), one(R), v
    end
    # p * u = q * v + r
    p, q, r = pseudo_divrem(u, v, algo)
    if iszero(q)
        not_divided_error(u, v)
    end
    if iszero(r)
        # Shortcut, does not change the output
        return zero(R), one(R), v
    end
    # TODO
    #if !algo.primitive_rem
    #    r = primitive_part(r, algo)::R
    #end
    # a * v + b * r = g
    # a * v + b * (p * u - q * v) = g
    # b * p * u + (a - b * q) * v = g
    a, b, g = primitive_univariate_gcdx(v, r, algo)
    return p * b, (a - b * q), g
end


function primitive_univariate_gcd(p::APL, q::APL, ::SubresultantAlgorithm)
    error("Not implemented yet")
end

"""
    univariate_gcd(p1::AbstractPolynomialLike, p2::AbstractPolynomialLike, algo::AbstractUnivariateGCDAlgorithm)

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
function univariate_gcd(p1::APL{S}, p2::APL{T}, algo::AbstractUnivariateGCDAlgorithm) where {S,T}
    return univariate_gcd(_field_absorb(algebraic_structure(S), algebraic_structure(T)), p1, p2, algo)
end
function univariate_gcd(::UFD, p1::APL, p2::APL, algo::AbstractUnivariateGCDAlgorithm)
    f1, g1 = primitive_part_content(p1, algo)
    f2, g2 = primitive_part_content(p2, algo)
    pp = primitive_univariate_gcd(f1, f2, algo)
    gg = _gcd(g1, g2, algo)#::MA.promote_operation(gcd, typeof(g1), typeof(g2))
    # Multiply each coefficient by the gcd of the contents.
    return mapcoefficientsnz(Base.Fix1(*, gg), pp)
end

function univariate_gcd(::Field, p1::APL, p2::APL, algo::AbstractUnivariateGCDAlgorithm)
    return primitive_univariate_gcd(p1, p2, algo)
end

function univariate_gcdx(p1::APL{S}, p2::APL{T}, algo::AbstractUnivariateGCDAlgorithm) where {S,T}
    return univariate_gcdx(_field_absorb(algebraic_structure(S), algebraic_structure(T)), p1, p2, algo)
end
function univariate_gcdx(::UFD, p1::APL, p2::APL, algo::AbstractUnivariateGCDAlgorithm)
    f1, g1 = primitive_part_content(p1, algo)
    f2, g2 = primitive_part_content(p2, algo)
    a, b, pp = primitive_univariate_gcdx(f1, f2, algo)
    gg = _gcd(g1, g2, algo)#::MA.promote_operation(gcd, typeof(g1), typeof(g2))
    # Multiply each coefficient by the gcd of the contents.
    return g2 * a, g1 * b, g1 * g2 * mapcoefficientsnz(Base.Fix1(*, gg), pp)
end
function univariate_gcdx(::Field, p1::APL, p2::APL, algo::AbstractUnivariateGCDAlgorithm)
    return primitive_univariate_gcdx(p1, p2, algo)
end

_gcd(a::APL, b::APL, algo) = gcd(a, b, algo)
_gcd(a, b::APL, algo) = gcd(a, b, algo)
_gcd(a::APL, b, algo) = gcd(a, b, algo)
_gcd(a, b, algo) = gcd(a, b)

_simplifier(a::APL, b::APL, algo) = gcd(a, b, algo)
_simplifier(a, b, algo) = _gcd(a, b, algo)
# Before Julia v1.4, it is not defined.
# After Julia v1.4, it is defined as `gcd of num / lcm of den`.
# We prefer `gcd of den`, otherwise,
# `1/a0 + 1/a1 x + 1/a2 x^2 + ... + 1/an x^n`
# will be transformed into
# `a1*a2*...*an + a0*a2*...*an x + ...`
# which makes the size of the `BigInt`s grow significantly which slows things down.
_simplifier(a::Rational, b::Rational, algo) = gcd(a.num, b.num) // gcd(a.den, b.den)

"""
    content(poly::AbstractPolynomialLike{T}, algo::AbstractUnivariateGCDAlgorithm) where {T}

Return the *content* of the polynomial `poly` over a unique factorization
domain `S` as defined in [Knu14, (3) p. 423].
That is, return the `gcd` of the coefficients of `poly`.
See also [`primitive_part_content`](@ref).

[Knu14] Knuth, D.E., 2014.
*Art of computer programming, volume 2: Seminumerical algorithms.*
Addison-Wesley Professional. Third edition.
"""
function content(poly::APL{T}, algo::AbstractUnivariateGCDAlgorithm) where {T}
    P = MA.promote_operation(gcd, T, T)
    # This is tricky to infer a `content` calls `gcd` which calls `content`, etc...
    # To help Julia break the loop, we annotate the result here.
    return reduce(
        (a, b) -> _simplifier(a, b, algo),
        coefficients(poly),
        init = zero(P),
    )::P
end
function content(poly::APL{T}, algo::AbstractUnivariateGCDAlgorithm) where {T<:AbstractFloat}
    return one(T)
end

"""
    primitive_part_content(poly::AbstractPolynomialLike{T}, algo::AbstractUnivariateGCDAlgorithm) where {T}

Return the *primitive part* of the polynomial `poly` over a unique
factorization domain `S` as defined in [Knu14, (3) p. 423].
That is, return the exact division of `poly` by its [`content`](@ref).
If the content is also needed, call [`primitive_part_content`](@ref)
instead.

[Knu14] Knuth, D.E., 2014.
*Art of computer programming, volume 2: Seminumerical algorithms.*
Addison-Wesley Professional. Third edition.
"""
primitive_part(p::APL, algo::AbstractUnivariateGCDAlgorithm) = primitive_part_content(p, algo)[1]
primitive_part(p::APL{<:AbstractFloat}, algo::AbstractUnivariateGCDAlgorithm) = p

"""
    primitive_part_content(poly::AbstractPolynomialLike{T}, algo::AbstractUnivariateGCDAlgorithm) where {T}

Return the *primitive part* and *content* of the polynomial `poly` over a unique
factorization domain `S` as defined in [Knu14, (3) p. 423]. This is more
efficient to call this function rather than calling [`primitive_part`](@ref) and
[`content`](@ref) separately since computing the primitive part requires
computing the content first and this function avoid computing the content twice.

[Knu14] Knuth, D.E., 2014.
*Art of computer programming, volume 2: Seminumerical algorithms.*
Addison-Wesley Professional. Third edition.
"""
function primitive_part_content(p, algo::AbstractUnivariateGCDAlgorithm)
    g = content(p, algo)
    return mapcoefficientsnz(Base.Fix2(_div, g), p), g
end

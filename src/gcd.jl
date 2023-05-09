export GeneralizedEuclideanAlgorithm

_copy(p, ::MA.IsMutable) = p
# `Base.copy` does not copy anything for `BigInt` so we need `MA.copy_if_mutable`
_copy(p, ::MA.IsNotMutable) = MA.copy_if_mutable(p)

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
function Base.gcd(α, p::APL, algo::AbstractUnivariateGCDAlgorithm=GeneralizedEuclideanAlgorithm(), mα::MA.MutableTrait=MA.IsNotMutable(), mp::MA.MutableTrait=MA.IsNotMutable())
    return _coefficient_gcd(α, content(p, algo, mp))
end
function Base.gcd(p::APL, α, algo::AbstractUnivariateGCDAlgorithm=GeneralizedEuclideanAlgorithm(), mp::MA.MutableTrait=MA.IsNotMutable(), mα::MA.MutableTrait=MA.IsNotMutable())
    return _coefficient_gcd(content(p, algo, mp), α)
end

#function MA.promote_operation(::typeof(gcd), P::Type{<:Number}, Q::Type{<:Number})
#    return typeof(gcd(one(P), one(Q)))
#end
function MA.promote_operation(::typeof(gcd), P::Type{<:APL}, Q::Type{<:APL}, A::Type=GeneralizedEuclideanAlgorithm)
    return MA.promote_operation(rem_or_pseudo_rem, P, Q, A)
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
function Base.gcd(p1::APL{T}, p2::APL{S}, algo::AbstractUnivariateGCDAlgorithm=GeneralizedEuclideanAlgorithm(), m1::MA.MutableTrait=MA.IsNotMutable(), m2::MA.MutableTrait=MA.IsNotMutable()) where {T, S}
    # If one of these is zero, `shift` should be infinite
    # for this method to work so we exclude these cases.
    if isapproxzero(p1)
        return convert(MA.promote_operation(gcd, typeof(p1), typeof(p2)), _copy(p2, m2))
    end
    if isapproxzero(p2)
        return convert(MA.promote_operation(gcd, typeof(p1), typeof(p2)), _copy(p1, m1))
    end
    shift1, defl1 = deflation(p1)
    shift2, defl2 = deflation(p2)
    shift = gcd(shift1, shift2)
    defl = map_exponents(gcd, defl1, defl2)
    # We factor out `x.^shift1` from `p1` and
    # `x.^shift2` from `p2`. The `gcd` of these
    # monomials is `x.^shift`.
    # Then, we subsitute `y[i] = x[i]^defl[i]`.
    q1 = deflate(p1, shift1, defl)
    q2 = deflate(p2, shift2, defl)
    g = deflated_gcd(q1, q2, algo, m1, m2)
    return inflate(g, shift, defl)::MA.promote_operation(gcd, typeof(p1), typeof(p2))
end

function Base.gcd(t1::AbstractTermLike{T}, t2::AbstractTermLike{S}, algo::AbstractUnivariateGCDAlgorithm=GeneralizedEuclideanAlgorithm(), m1::MA.MutableTrait=MA.IsNotMutable(), m2::MA.MutableTrait=MA.IsNotMutable()) where {T, S}
    return term(_coefficient_gcd(coefficient(t1), coefficient(t2)), gcd(monomial(t1), monomial(t2)))
end

function shift_deflation(p::AbstractPolynomialLike, v::AbstractVariable)
    shift = -1
    defl = 0
    for mono in monomials(p)
        exp = degree(mono, v)
        if shift == -1
            shift = exp
        elseif exp < shift
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
            defl = gcd(defl, shift - exp)
            shift = exp
        else
            defl = gcd(defl, exp - shift)
        end
    end
    return shift, defl
end

# Inspired from to `AbstractAlgebra.deflation`
function deflation(p::AbstractPolynomialLike)
    if iszero(p)
        return constant_monomial(p), constant_monomial(p)
    end
    shift_defl = shift_deflation.(p, variables(p))
    shift = getindex.(shift_defl, 1)
    defl = getindex.(shift_defl, 2)
    s = prod(variables(p).^shift; init = constant_monomial(p))::monomial_type(p)
    d = prod(variables(p).^defl; init = constant_monomial(p))::monomial_type(p)
    return s, d
end

function _zero_to_one_exp(defl::AbstractMonomial)
    # TODO Make it faster by calling something like `map_exponents`.
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
    return inflate(convert(polynomial_type(shift, typeof(α)), α), shift, defl)
end
function inflate(p::AbstractPolynomialLike, shift, defl)
    if isconstant(shift) && all(d -> isone(d) || iszero(d), exponents(defl))
        return p
    end
    q = MA.operate(inflate, p, shift, _zero_to_one_exp(defl))
    return q
end

function MA.operate(::typeof(deflate), mono::AbstractMonomial, shift, defl)
    mutable_mono = map_exponents(-, mono, shift)
    return map_exponents!(div, mutable_mono, defl)
end
function MA.operate(::typeof(inflate), mono::AbstractMonomial, shift, defl)
    mutable_mono = map_exponents(*, mono, defl)
    return map_exponents!(+, mutable_mono, shift)
end

# Inspired from to `AbstractAlgebra.deflate`
function MA.operate(op::Union{typeof(deflate), typeof(inflate)}, p::AbstractPolynomialLike, shift, defl)
    return polynomial(map(terms(p)) do t
        return term(coefficient(t), MA.operate(op, monomial(t), shift, defl))
    end)
end

function deflated_gcd(p1::APL{T}, p2::APL{S}, algo, m1::MA.MutableTrait, m2::MA.MutableTrait) where {T, S}
    i1, i2, num_common = _extracted_variable(p1, p2)
    if iszero(i1)
        if iszero(i2)
            return univariate_gcd(p1, p2, algo, m1, m2)
        else
            if isapproxzero(p1)
                return convert(MA.promote_operation(gcd, typeof(p1), typeof(p2)), _copy(p2, m2))
            end
            v2 = variables(p2)[i2]
            q2 = isolate_variable(p2, v2, m2)
            g = content(q2, algo, MA.IsMutable())
            return gcd(p1, g, algo, m1, MA.IsMutable())
        end
    else
        if iszero(i2)
            if isapproxzero(p2)
                return convert(MA.promote_operation(gcd, typeof(p1), typeof(p2)), _copy(p1, m1))
            end
            v1 = variables(p1)[i1]
            q1 = isolate_variable(p1, v1, m1)
            g = content(q1, algo, MA.IsMutable())
            return gcd(g, p2, algo, MA.IsMutable(), m2)
        else
            if num_common > 1
                v1 = variables(p1)[i1]
                @assert v1 == variables(p2)[i2]
                return multivariate_gcd(p1, p2, v1, algo, m1, m2)
            else
                return univariate_gcd(p1, p2, algo, m1, m2)
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

function multivariate_gcd(p1::APL, p2::APL, var, algo, m1::MA.MutableTrait, m2::MA.MutableTrait)
    q1 = isolate_variable(p1, var, m1)
    q2 = isolate_variable(p2, var, m2)
    q = univariate_gcd(q1, q2, algo, MA.IsMutable(), MA.IsMutable())
    P = MA.promote_operation(gcd, typeof(p1), typeof(p2))
    return flatten_variable!(term_type(P), q)::P
end

_vector(t::AbstractVector) = collect(t)
_vector(t::Vector) = t

"""
    isolate_variable(poly::APL, var::AbstractVariable, mutability::MA.MutableTrait)

Returns a polynomial with variable `var`. The other variables of `poly` are moved as coefficients.

The output can be mutated without affecting `poly` if `mutability` is
`MA.IsNotMutable`.
"""
function isolate_variable(poly::APL, var::AbstractVariable, mutability::MA.MutableTrait)
    old_terms = sort!(_vector(terms(_copy(poly, mutability))), by = Base.Fix2(degree, var))
    U = MA.promote_operation(substitute, Subs, typeof(poly), Pair{typeof(var),Int})
    T = term_type(var, U)
    new_terms = T[]
    i = firstindex(old_terms)
    while i <= lastindex(old_terms)
        j = i + 1
        d = degree(old_terms[i], var)
        while j <= lastindex(old_terms)
            if degree(old_terms[j], var) != d
                break
            end
            j += 1
        end
        coef = _polynomial([subs(old_terms[k], (var,) => (1,)) for k in i:(j - 1)], SortedUniqState(), mutability)
        push!(new_terms, term(coef, var^d))
        i = j
    end
    return polynomial!(new_terms, SortedUniqState())
end

function flatten_variable!(::Type{TT}, poly::APL) where {TT<:AbstractTerm}
    ts = TT[]
    for t in terms(poly)
        m = monomial(t)
        for _t in terms(coefficient(t))
            push!(ts, _t * m)
        end
    end
    return polynomial!(ts, UniqState())
end

_polynomial(ts, state, ::MA.IsNotMutable) = polynomial(ts, state)
_polynomial(ts, state, ::MA.IsMutable) = polynomial!(ts, state)

"""
    primitive_univariate_gcd!(p::APL, q::APL, algo::AbstractUnivariateGCDAlgorithm)

Returns the `gcd` of primitive polynomials `p` and `q` using algorithm `algo`
which is a subtype of [`AbstractUnivariateGCDAlgorithm`](@ref).
The function might modify `p` or `q`.
"""
function primitive_univariate_gcd! end

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

# If `p` and `q` do not have the same type then the local variables `p` and `q`
# won't be type stable so we create `u` and `v`.
function primitive_univariate_gcd!(p::APL, q::APL, algo::GeneralizedEuclideanAlgorithm)
    if maxdegree(p) < maxdegree(q)
        return primitive_univariate_gcd!(q, p, algo)
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
            return MA.operate!(one, u)
        end

        d_before = degree(leading_monomial(u))
        r = MA.operate!(rem_or_pseudo_rem, u, v, algo)
        d_after = degree(leading_monomial(r))
        if d_after == d_before
            not_divided_error(u, v)
        end

        if !algo.primitive_rem
            r = primitive_part(r, algo, MA.IsMutable())::R
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


function primitive_univariate_gcd!(p::APL, q::APL, ::SubresultantAlgorithm)
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
   [`primitive_univariate_gcd!`](@ref) for primitive parts.
3. Return the product of these two `gcd`; see
   [Knu14, Algorithm E: E4, p. 427] or [Knu14, Algorithm C: C4, p. 429].

[Knu14] Knuth, D.E., 2014.
*Art of computer programming, volume 2: Seminumerical algorithms.*
Addison-Wesley Professional. Third edition.
"""
function univariate_gcd(p1::APL{S}, p2::APL{T}, algo::AbstractUnivariateGCDAlgorithm, m1::MA.MutableTrait, m2::MA.MutableTrait) where {S,T}
    return univariate_gcd(_field_absorb(algebraic_structure(S), algebraic_structure(T)), p1, p2, algo, m1, m2)
end
function univariate_gcd(::UFD, p1::APL, p2::APL, algo::AbstractUnivariateGCDAlgorithm, m1::MA.MutableTrait, m2::MA.MutableTrait)
    f1, g1 = primitive_part_content(p1, algo, m1)
    f2, g2 = primitive_part_content(p2, algo, m2)
    pp = primitive_univariate_gcd!(f1, f2, algo)
    gg = _gcd(g1, g2, algo, MA.IsMutable(), MA.IsMutable())#::MA.promote_operation(gcd, typeof(g1), typeof(g2))
    # Multiply each coefficient by the gcd of the contents.
    if !isone(gg)
        MA.operate!(right_constant_mult, pp, gg)
    end
    return pp
end

function univariate_gcd(::Field, p1::APL, p2::APL, algo::AbstractUnivariateGCDAlgorithm, m1::MA.MutableTrait, m2::MA.MutableTrait)
    return primitive_univariate_gcd!(_copy(p1, m1), _copy(p2, m2), algo)
end

function univariate_gcdx(p1::APL{S}, p2::APL{T}, algo::AbstractUnivariateGCDAlgorithm) where {S,T}
    return univariate_gcdx(_field_absorb(algebraic_structure(S), algebraic_structure(T)), p1, p2, algo)
end
function univariate_gcdx(::UFD, p1::APL, p2::APL, algo::AbstractUnivariateGCDAlgorithm)
    f1, g1 = primitive_part_content(p1, algo, MA.IsNotMutable())
    f2, g2 = primitive_part_content(p2, algo, MA.IsNotMutable())
    a, b, pp = primitive_univariate_gcdx(f1, f2, algo)
    gg = _gcd(g1, g2, algo, MA.IsMutable(), MA.IsMutable())#::MA.promote_operation(gcd, typeof(g1), typeof(g2))
    # Multiply each coefficient by the gcd of the contents.
    return g2 * a, g1 * b, g1 * g2 * map_coefficients(Base.Fix1(*, gg), pp, nonzero = true)
end
function univariate_gcdx(::Field, p1::APL, p2::APL, algo::AbstractUnivariateGCDAlgorithm)
    return primitive_univariate_gcdx(p1, p2, algo)
end

_gcd(a::APL, b::APL, algo, ma, mb) = gcd(a, b, algo, ma, mb)
_gcd(a, b::APL, algo, ma, mb) = gcd(a, b, algo, ma, mb)
_gcd(a::APL, b, algo, ma, mb) = gcd(a, b, algo, ma, mb)
_gcd(a, b, algo, ma, mb) = gcd(a, b)

_simplifier(a::APL, b::APL, algo, ma, mb) = gcd(a, b, algo, ma, mb)
_simplifier(a, b, algo, ma, mb) = _gcd(a, b, algo, ma, mb)
# Before Julia v1.4, it is not defined.
# After Julia v1.4, it is defined as `gcd of num / lcm of den`.
# We prefer `gcd of den`, otherwise,
# `1/a0 + 1/a1 x + 1/a2 x^2 + ... + 1/an x^n`
# will be transformed into
# `a1*a2*...*an + a0*a2*...*an x + ...`
# which makes the size of the `BigInt`s grow significantly which slows things down.
_simplifier(a::Rational, b::Rational, algo, ma, mb) = gcd(a.num, b.num) // gcd(a.den, b.den)

# Largely inspired from from `YingboMa/SIMDPolynomials.jl`.
function termwise_content(p::APL, algo, mutability::MA.MutableTrait)
    ts = terms(p)
    length(ts) == 1 && return _copy(first(ts), mutability)
    g = gcd(ts[1], ts[2], algo, mutability, mutability)
    isone(g) || for i in 3:length(ts)
        g = gcd(g, ts[i], algo, MA.IsMutable(), mutability)
        isone(g) && break
    end
    return g
end

"""
    content(poly::AbstractPolynomialLike{T}, algo::AbstractUnivariateGCDAlgorithm, mutability::MA.MutableTrait) where {T}

Return the *content* of the polynomial `poly` over a unique factorization
domain `S` as defined in [Knu14, (3) p. 423].
That is, return the `gcd` of the coefficients of `poly`.
See also [`primitive_part_content`](@ref).

The output can be mutated without affecting `poly` if `mutability` is
`MA.IsNotMutable`.

[Knu14] Knuth, D.E., 2014.
*Art of computer programming, volume 2: Seminumerical algorithms.*
Addison-Wesley Professional. Third edition.
"""
function content(poly::APL{T}, algo::AbstractUnivariateGCDAlgorithm, mutability::MA.MutableTrait) where {T}
    P = MA.promote_operation(gcd, T, T)
    coefs = coefficients(poly)
    if isempty(coefs)
        return zero(P)
    end
    if length(coefs) == 1
        return convert(P, _copy(first(coefs), mutability))
    end
    # Largely inspired from from `YingboMa/SIMDPolynomials.jl`.
    if T <: APL
        for i in eachindex(coefs)
            if nterms(coefs[i]) == 1
                g = _gcd(termwise_content(coefs[1], algo, mutability), termwise_content(coefs[2], algo, mutability), algo, MA.IsMutable(), MA.IsMutable())
                isone(g) || for i in 3:length(coefs)
                    g = _gcd(g, termwise_content(coefs[i], algo, mutability), algo, MA.IsMutable(), MA.IsMutable())
                    isone(g) && break
                end
                return convert(P, g)
            end
        end
    end
    # This is tricky to infer a `content` calls `gcd` which calls `content`, etc...
    # To help Julia break the loop, we annotate the result here.
    g = _gcd(coefs[1], coefs[2], algo, mutability, mutability)::P
    isone(g) || for i in 3:length(coefs)
        g = _simplifier(g, coefs[i], algo, MA.IsMutable(), mutability)::P
        isone(g) && break
    end
    return g::P
end
function content(::APL{T}, ::AbstractUnivariateGCDAlgorithm, ::MA.MutableTrait) where {T<:AbstractFloat}
    return one(T)
end

"""
    primitive_part(poly::AbstractPolynomialLike{T}, algo::AbstractUnivariateGCDAlgorithm) where {T}

Return the *primitive part* of the polynomial `poly` over a unique
factorization domain `S` as defined in [Knu14, (3) p. 423].
That is, return the exact division of `poly` by its [`content`](@ref).
If the content is also needed, call [`primitive_part_content`](@ref)
instead.

[Knu14] Knuth, D.E., 2014.
*Art of computer programming, volume 2: Seminumerical algorithms.*
Addison-Wesley Professional. Third edition.
"""
primitive_part(p::APL, algo::AbstractUnivariateGCDAlgorithm, mutability::MA.MutableTrait) = primitive_part_content(p, algo, mutability)[1]
primitive_part(p::APL{<:AbstractFloat}, ::AbstractUnivariateGCDAlgorithm, ::MA.MutableTrait) = p

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
function primitive_part_content(p, algo::AbstractUnivariateGCDAlgorithm, mutability::MA.MutableTrait)
    g = content(p, algo, MA.IsNotMutable())
    return right_constant_div_multiple(p, g, mutability), g
end

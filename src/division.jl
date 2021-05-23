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
_div(a::Union{Rational, AbstractFloat}, b) = a / b
_div(a, b) = div(a, b)
_div(m1::AbstractMonomialLike, m2::AbstractMonomialLike) = mapexponents(-, m1, m2)
function _div(t::AbstractTerm, m::AbstractMonomial)
    term(coefficient(t), _div(monomial(t), m))
end
function _div(t1::AbstractTermLike, t2::AbstractTermLike)
    term(_div(coefficient(t1), coefficient(t2)), _div(monomial(t1), monomial(t2)))
end
function _div(f::APL, g::APL)
    lt = leadingterm(g)
    rf = MA.copy_if_mutable(f)
    rg = removeleadingterm(g)
    q = zero(rf)
    while !iszero(rf)
        ltf = leadingterm(rf)
        qt = _div(ltf, lt)
        q = MA.add!(q, qt)
        rf = MA.operate!(removeleadingterm, rf)
        rf = MA.operate!(MA.sub_mul, rf, qt, rg)
    end
    return q
end

Base.div(f::APL, g::Union{APL, AbstractVector{<:APL}}; kwargs...) = divrem(f, g; kwargs...)[1]
Base.rem(f::APL, g::Union{APL, AbstractVector{<:APL}}; kwargs...) = divrem(f, g; kwargs...)[2]

# FIXME What should we do for `Rational` ?
function ring_rem(f::APL, g::APL{<:Union{Rational, AbstractFloat}})
    return true, rem(f, g)
end
function ring_rem(f::APL, g::APL)
    ltg = leadingterm(g)
    ltf = leadingterm(f)
    if !divides(monomial(ltg), ltf)
        return false, f
    end
    new_f = constantterm(coefficient(ltg), f) * removeleadingterm(f)
    new_g = term(coefficient(ltf), _div(monomial(ltf), monomial(ltg))) * removeleadingterm(g)
    return true, constantterm(coefficient(ltg), f) * removeleadingterm(f) -
        term(coefficient(ltf), _div(monomial(ltf), monomial(ltg))) * removeleadingterm(g)
end
function MA.promote_operation(
    ::typeof(ring_rem),
    ::Type{P},
    ::Type{Q},
) where {T,S<:Union{Rational, AbstractFloat},P<:APL{T},Q<:APL{S}}
    return MA.promote_operation(rem, P, Q)
end
function MA.promote_operation(::typeof(ring_rem), ::Type{P}, ::Type{Q}) where {T,S,P<:APL{T},Q<:APL{S}}
    U1 = MA.promote_operation(*, S, T)
    U2 = MA.promote_operation(*, T, S)
    # `promote_type(P, Q)` is needed for TypedPolynomials in case they use different variables
    return polynomialtype(promote_type(P, Q), MA.promote_operation(-, U1, U2))
end

function MA.promote_operation(::Union{typeof(div), typeof(rem)}, ::Type{P}, g::Type{Q}) where {T,S,P<:APL{T},Q<:APL{S}}
    U = MA.promote_operation(/, T, S)
    # `promote_type(P, Q)` is needed for TypedPolynomials in case they use different variables
    return polynomialtype(promote_type(P, Q), MA.promote_operation(-, U, U))
end
function Base.divrem(f::APL{T}, g::APL{S}; kwargs...) where {T, S}
    rf = convert(MA.promote_operation(div, typeof(f), typeof(g)), MA.copy_if_mutable(f))
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
    rf = convert(MA.promote_operation(div, typeof(f), eltype(g)), MA.copy_if_mutable(f))
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

_simplifier(a, b) = gcd(a, b)
# Before Julia v1.4, it is not defined.
# After Julia v1.4, it is defined as `gcd of num / lcm of den`.
# We prefer `gcd of den`, otherwise,
# `1/a0 + 1/a1 x + 1/a2 x^2 + ... + 1/an x^n`
# will be transformed into
# `a1*a2*...*an + a0*a2*...*an x + ...`
# which makes the size of the `BigInt`s grow significantly which slows things down.
_simplifier(a::Rational, b::Rational) = gcd(a.num, b.num) // gcd(a.den, b.den)

function coefficients_gcd(poly::APL{T}) where {T}
    P = MA.promote_operation(gcd, T, T)
    # This is tricky to infer a `coefficients_gcd` calls `gcd` which calls `coefficients_gcd`, etc...
    # To help Julia break the loop, we annotate the result here.
    return reduce(
        _simplifier,
        coefficients(poly),
        init = zero(P),
    )::P
end

#function MA.promote_operation(::typeof(gcd), P::Type{<:Number}, Q::Type{<:Number})
#    return typeof(gcd(one(P), one(Q)))
#end
function MA.promote_operation(::typeof(gcd), P::Type{<:APL}, Q::Type{<:APL})
    return MA.promote_operation(ring_rem, P, Q)
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
    s = prod(variables(p).^shift)
    d = prod(variables(p).^defl)
    return s, d
end

function deflate(p::AbstractPolynomialLike, shift, defl)
    if isconstant(shift) && all(d -> isone(d) || iszero(d), exponents(defl))
        return p
    end
    q = MA.operate(deflate, p, shift, defl)
    return q
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
    mapexponents!(div, mutable_mono, defl)
    return mutable_mono
end
function MA.operate(::typeof(inflate), mono::AbstractMonomial, shift, defl)
    mutable_mono = mapexponents(*, mono, defl)
    mapexponents!(+, mutable_mono, shift)
    return mutable_mono
end

# Inspired from to `AbstractAlgebra.deflate`
function MA.operate(op::Union{typeof(deflate), typeof(inflate)}, p::AbstractPolynomialLike, shift, defl)
    defl = prod(variables(defl).^map(d -> iszero(d) ? one(d) : d, exponents(defl)))
    return polynomial(map(terms(p)) do t
        return term(coefficient(t), MA.operate(op, monomial(t), shift, defl))
    end)
end

"""
    function gcd(p1::AbstractPolynomialLike{T}, p2::AbstractPolynomialLike{S}) where {T, S}

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
    return inflate(g, shift, defl)
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
            g = coefficients_gcd(q2)
            return gcd(p1, g)
        end
    else
        if iszero(i2)
            if isapproxzero(p2)
                return convert(MA.promote_operation(gcd, typeof(p1), typeof(p2)), p1)
            end
            v1 = variables(p1)[i1]
            q1 = isolate_variable(p1, v1)
            g = coefficients_gcd(q1)
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

_div_gcd_coefficients(p::APL{<:AbstractFloat}) = p
_div_gcd_coefficients(p::APL) = div_gcd_coefficients(p)[1]

function _gcd_relatively_prime_coefficients(p::APL, q::APL)
    if isapproxzero(q)
        convert(MA.promote_operation(gcd, typeof(p), typeof(q)), p)
    elseif isapproxzero(p)
        convert(MA.promote_operation(gcd, typeof(p), typeof(q)), q)
    else
        divided, r = ring_rem(p, q)
        o = q
        if !divided
            divided, r = ring_rem(q, p)
            o = p
            # Since `p` and `q` are univariate, at least one divides the other
            @assert divided
        end
        return _gcd_relatively_prime_coefficients(o, _div_gcd_coefficients(r))
    end
end
function univariate_gcd(p1::APL, p2::APL{<:AbstractFloat})
    return _gcd_relatively_prime_coefficients(p1, p2)
end
function univariate_gcd(p1::APL, p2::APL)
    f1, g1 = div_gcd_coefficients(p1)
    f2, g2 = div_gcd_coefficients(p2)
    pp = _gcd_relatively_prime_coefficients(f1, f2)
    gg = gcd(g1, g2)#::MA.promote_operation(gcd, typeof(g1), typeof(g2))
    return mapcoefficientsnz(Base.Fix1(*, gg), pp)
end
function div_gcd_coefficients(p)
    g = coefficients_gcd(p)
    return mapcoefficientsnz(Base.Fix2(_div, g), p), g
end
function multivariate_gcd(p1::APL, p2::APL, var)
    q = univariate_gcd(isolate_variable(p1, var), isolate_variable(p2, var))
    return sum(coefficient(t) * monomial(t) for t in terms(q))::MA.promote_operation(gcd, typeof(p1), typeof(p2))
end

Base.lcm(p::APL, q::APL) = p * div(q, gcd(p, q))

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

function Base.gcd(p::AbstractPolynomialLike{T}, q::AbstractPolynomialLike{S}) where {T, S}
    if isapproxzero(q)
        convert(polynomialtype(p, Base.promote_op(proddiff, T, S)), p)
    else
        gcd(q, rem(p, q))
    end
end
Base.lcm(p::AbstractPolynomialLike, q::AbstractPolynomialLike) = p * div(q, gcd(p, q))

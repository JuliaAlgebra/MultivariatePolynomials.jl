"""
    divides(t1::AbstractTermLike, t2::AbstractTermLike)

Returns whether the monomial of t1 divides the monomial of t2.

# Examples

Calling `divides(3x*y, 2x^2*y)` should return true but calling `divides(x*y, y)` should return false.
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
function _div(t1::AbstractTerm, t2::AbstractTerm)
    (coefficient(t1) / coefficient(t2)) * _div(monomial(t1), monomial(t2))
end

proddiff(x, y) = x*y - x*y

function Base.divrem(f::APL{T}, g::APL{S}) where {T, S}
    rf = changecoefficienttype(f, Base.promote_op(proddiff, T, S))
    q = r = zero(f - g / 2)
    lt = leadingterm(g)
    rg = removeleadingterm(g)
    lm = monomial(lt)
    while !iszero(rf)
        ltf = leadingterm(rf)
        if divides(lm, ltf)
            qt = _div(ltf, lt)
            q += qt
            rf = removeleadingterm(rf) - qt * rg
        elseif lm > monomial(ltf)
            # Since the monomials are sorted in decreasing order,
            # lm is larger than all of them hence it cannot divide any of them
            r += rf
            break
        else
            r += ltf
            rf = removeleadingterm(rf)
        end
    end
    q, r
end
Base.div(f::APL, g::APL) = divrem(f, g)[1]
Base.rem(f::APL, g::APL) = divrem(f, g)[2]

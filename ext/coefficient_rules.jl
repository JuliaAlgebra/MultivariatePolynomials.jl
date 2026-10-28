# Normalize ChainRules cotangents before applying the algebraic pullback.
function _rrule(f, args...)
    y, pb = _coefficient_rrule(f, args...)
    function pullback(delta)
        d = CRC.unthunk(delta)
        return pb(d isa CRC.AbstractZero ? zero(y) : d)
    end
    return y, pullback
end

# The Hermitian coefficient pairing; reverse AD uses its real part.
# This is independent of MP's polynomial-valued dot product.
function coefficient_dot(p, q)
    T = MA.promote_operation(
        *,
        MA.promote_operation(conj, coefficient_type(p)),
        coefficient_type(q),
    )
    return sum(terms(p); init = zero(T)) do t
        return conj(coefficient(t)) * coefficient(q, monomial(t))
    end
end
coefficient_at(p::Number, m) = isone(m) ? p : zero(p)
coefficient_at(p, m) = coefficient(p, m)

function adjoint_multiply(d, q, y)
    d isa CRC.AbstractZero && return d
    d = d isa Number ? d * one(y) : d
    T = promote_type(coefficient_type(d), coefficient_type(q))
    result = zero(T, parent(polynomial(d)))
    for t in terms(d), s in terms(q)
        if divides(monomial(s), monomial(t))
            result += term(
                coefficient(t) * conj(coefficient(s)),
                MP.div_multiple(monomial(t), monomial(s)),
            )
        end
    end
    return result
end

function _coefficient_rrule(::typeof(*), p::MP._APL, q::MP._APL)
    is_commutative(p) && is_commutative(q) || throw(
        ArgumentError("Polynomial AD currently requires commutative monomials"),
    )
    y = p * q
    pb(d) =
        (CRC.NoTangent(), adjoint_multiply(d, q, y), adjoint_multiply(d, p, y))
    return y, pb
end
for op in (:+, :-)
    @eval function _coefficient_rrule(::typeof($op), p::MP._APL, c::Number)
        return $op(p, c),
        d -> (CRC.NoTangent(), d, $op(coefficient_at(d, constant_monomial(p))))
    end
    @eval function _coefficient_rrule(::typeof($op), c::Number, p::MP._APL)
        return $op(c, p),
        d -> (CRC.NoTangent(), coefficient_at(d, constant_monomial(p)), $op(d))
    end
end
function _coefficient_rrule(::typeof(*), p::MP._APL, c::Number)
    y = p * c
    pb(d) = (
        CRC.NoTangent(),
        d * conj(c),
        coefficient_dot(p, d isa Number ? d * one(y) : d),
    )
    return y, pb
end
function _coefficient_rrule(::typeof(*), c::Number, p::MP._APL)
    y, pb = _coefficient_rrule(*, p, c)
    function pullback(d)
        _, dp, dc = pb(d)
        return CRC.NoTangent(), dc, dp
    end
    return y, pullback
end
function _coefficient_rrule(::typeof(/), p::MP._APL, c::Number)
    y = p / c
    pb(d) = (
        CRC.NoTangent(),
        d / conj(c),
        -coefficient_dot(y, d isa Number ? d * one(y) : d) / conj(c),
    )
    return y, pb
end
function _coefficient_rrule(::typeof(^), p::MP._APL, n::Integer)
    is_commutative(p) || throw(
        ArgumentError("Polynomial AD currently requires commutative monomials"),
    )
    y = p^n
    pb(d) = (
        CRC.NoTangent(),
        iszero(n) ? zero(y) : adjoint_multiply(d, n * p^(n - 1), y),
        CRC.NoTangent(),
    )
    return y, pb
end
function _coefficient_rrule(
    ::typeof(differentiate),
    p::MP._APL,
    x::AbstractVariable,
)
    is_commutative(p) || throw(
        ArgumentError("Polynomial AD currently requires commutative monomials"),
    )
    y = differentiate(p, x)
    pb(d) = (CRC.NoTangent(), x * differentiate(x * d, x), CRC.NoTangent())
    return y, pb
end
function _coefficient_rrule(
    ::typeof(coefficient),
    p::MP._APL,
    m::AbstractMonomialLike,
)
    return coefficient(p, m), d -> (CRC.NoTangent(), d * m, CRC.NoTangent())
end
function _coefficient_rrule(::typeof(coefficient), p::AbstractTerm)
    return coefficient(p), d -> (CRC.NoTangent(), d * monomial(p))
end
function _coefficient_rrule(::typeof(term), c::Number, m::AbstractMonomialLike)
    return term(c, m),
    d -> (CRC.NoTangent(), coefficient_at(d, m), CRC.NoTangent())
end
function _coefficient_rrule(
    ::typeof(polynomial),
    c::AbstractVector{<:Number},
    m::AbstractVector{<:AbstractMonomialLike},
)
    return polynomial(c, m),
    d -> (CRC.NoTangent(), [coefficient_at(d, x) for x in m], CRC.NoTangent())
end
function _coefficient_rrule(::typeof(polynomial), p::MP._APL)
    return polynomial(p), d -> (CRC.NoTangent(), d)
end

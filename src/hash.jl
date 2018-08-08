_hashpowers(u::UInt) = u
function _hashpowers(u::UInt, power::Tuple, powers...)
    if iszero(power[2])
        _hashpowers(u, powers...)
    else
        _hashpowers(hash(power, u), powers...)
    end
end
function Base.hash(m::AbstractMonomial, u::UInt)
    nnz = count(!iszero, exponents(m))
    if iszero(nnz)
        hash(1, u)
    elseif isone(nnz) && isone(degree(m))
        hash(variable(m), u)
    else
        _hashpowers(u, powers(m)...)
    end
end

function Base.hash(t::AbstractTerm, u::UInt)
    if iszero(t)
        hash(0, u)
    elseif coefficient(t) == 1
        hash(monomial(t), u)
    else
        hash(monomial(t), hash(coefficient(t), u))
    end
end

function Base.hash(p::AbstractPolynomial, u::UInt)
    if iszero(p)
        hash(0, u)
    else
        reduce((u, t) -> hash(t, u), terms(p); init=u)
    end
end

function Base.hash(p::Polynomial{B}, u::UInt) where {B}
    @assert is_commutative(variables(p))
    return hash(
        Tuple(
            (v, e) for (v, e) in zip(variables(p), exponents(p)) if !iszero(e)
        ),
        hash(B, u),
    )
end
function Base.hash(m::AbstractMonomial, u::UInt)
    @assert is_commutative(m)
    nnz = count(!iszero, exponents(m))
    if iszero(nnz)
        hash(1, u)
    elseif isone(nnz) && isone(degree(m))
        # We want the has of `m` to match the hash of `variable(m)`
        # so we ignore the exponent one.
        hash(variable(m), u)
    else
        reduce(
            (u, power) -> iszero(power[2]) ? u : hash(power, u),
            powers(m),
            init = u,
        )
    end
end

function Base.hash(t::AbstractTerm, u::UInt)
    if !iszero(t) && isconstant(t)
        return hash(coefficient(t), u)
    end
    return invoke(hash, Tuple{SA.Term,UInt}, t, u)
end

function Base.hash(p::AbstractPolynomial, u::UInt)
    if iszero(p)
        hash(0, u)
    else
        reduce((u, t) -> hash(t, u), terms(p); init = u)
    end
end

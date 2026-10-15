variables(b::_NoncommutativeBasis) = b.variables
nvariables(b::_NoncommutativeBasis) = length(variables(b))
monomial_type(::Type{<:_NoncommutativeBasis{M}}) where {M} = M
monomial_type(b::_NoncommutativeBasis) = monomial_type(typeof(b))
function Base.one(b::_NoncommutativeBasis)
    return monomial(variables(b), map(_ -> 0, variables(b)))
end
Base.isempty(::_NoncommutativeBasis) = false

function Base.:(==)(a::_NoncommutativeBasis, b::_NoncommutativeBasis)
    return variables(a) == variables(b)
end
Base.hash(b::_NoncommutativeBasis, h::UInt) = hash(variables(b), h)

function Base.in(m::AbstractMonomial, b::_NoncommutativeBasis)
    return all(v in variables(b) for (v, e) in powers(m) if !iszero(e))
end
Base.haskey(b::_NoncommutativeBasis, m::AbstractMonomial) = m in b
function Base.getindex(b::_NoncommutativeBasis{M}, m::M) where {M}
    m in b ||
        throw(ArgumentError("Monomial $m is outside the destination basis"))
    return m
end

Base.IteratorSize(::Type{<:_NoncommutativeBasis}) = Base.SizeUnknown()
function Base.iterate(b::_NoncommutativeBasis)
    return one(b), (1, monomials(variables(b), 1), 1)
end
function Base.iterate(b::_NoncommutativeBasis, (degree, monos, i))
    if i > length(monos)
        isempty(variables(b)) && return nothing
        degree += 1
        monos = monomials(variables(b), degree)
        i = 1
    end
    return monos[i], (degree, monos, i + 1)
end

SA.comparable(::_NoncommutativeBasis) = isless
_comparable_type(::Type{<:_NoncommutativeBasis}) = typeof(isless)
function keys_as_monomials(keys, ::_NoncommutativeBasis)
    return collect(keys)
end

function explicit_basis_covering(
    full::_NoncommutativeBasis,
    target::SubBasis{Monomial},
)
    @assert full == parent(target)
    return SA.SubBasis(full, target.keys)
end

function maxdegree_basis(basis::_NoncommutativeBasis, maxdeg::Int)
    return unsafe_basis(Monomial, monomials(variables(basis), 0:maxdeg))
end

function polynomial_type(
    ::Type{<:_NoncommutativeBasis{M}},
    ::Type{T},
) where {M,T}
    return polynomial_type(M, T)
end
function polynomial_type(b::_NoncommutativeBasis, ::Type{T}) where {T}
    return polynomial_type(typeof(b), T)
end

algebra(b::_NoncommutativeBasis) = SA.StarAlgebra(b, b)
function MA.promote_operation(
    ::typeof(algebra),
    ::Type{B},
) where {M,B<:_NoncommutativeBasis{M}}
    return SA.StarAlgebra{B,M,SA.DiracMStructure{M,M,B,typeof(*)}}
end

function SA.promote_bases_with_maps(a::B, b::B) where {B<:_NoncommutativeBasis}
    a == b && return (a, nothing), (b, nothing)
    vars = SA.merge_sorted(
        variables(a),
        variables(b);
        lt = isless,
        combine = SA.first_of,
        filter = _ -> true,
        rev = true,
    )
    basis = B(vars)
    # The parent may grow, but keys already carry their complete factor order.
    return (basis, identity), (basis, identity)
end
function SA.promote_object(::_NoncommutativeBasis, m::SA.DiracMStructure, _)
    return SA.basis(m)
end

# From MultivariateBases/src/bases.jl

const FullBasis{B,V,E} = SA.MappedBasis{Polynomial{B,V,E}}
const SubBasis{B,V,E} = SA.SubBasis{Polynomial{B,V,E}}
const MonomialIndexedBasis{B,V,E} = Union{SubBasis{B,V,E},FullBasis{B,V,E}}

exponents_index(::FullBasis, exp) = exp
exponents_index(b::SubBasis, exp) = SA.key_index(b, exp)

function ordering(::Type{MonomialIndexedBasis{B,V,E}}) where {B,V,E}
    return ordering(E)
end
ordering(b::MonomialIndexedBasis) = ordering(typeof(b))

function FullBasis(vars::Variables{B,V}) where {B,V}
    if !is_commutative(V)
        B === Monomial ||
            error("Noncommutative variables require the Monomial basis")
        return _NoncommutativeBasis{monomial_type(V),V}(
            unique!(sort(vars.variables; rev = true)),
        )
    end
    O = ordering(vars.variables)
    exps = ExponentsIterator{O}(constant_monomial_exponents(vars))
    return SA.MappedBasis{Polynomial{B,V,eltype(exps)}}(exps, vars, exponents)
end

function FullBasis{B}(vars) where {B}
    return FullBasis(Variables{B,typeof(vars)}(vars))
end

function FullBasis{B}(p::_APL) where {B}
    return FullBasis{B}(variables(p))
end
# Also accept Polynomial basis elements (monomials in the new world)
function FullBasis{B}(p::Polynomial) where {B}
    return FullBasis{B}(variables(p))
end

function SA.star(b::MonomialIndexedBasis{B,V,E}, exp::E) where {B,V,E}
    return b[SA.star(b[exp])]
end

constant_monomial_exponents(b::FullBasis) = constant_monomial_exponents(b.map)

monomial_type(::Type{<:FullBasis{B,V}}) where {B,V} = monomial_type(V)
function polynomial_type(basis::FullBasis, ::Type{T}) where {T}
    return polynomial_type(typeof(basis), T)
end

nvariables(basis::FullBasis) = nvariables(basis.map)
nvariables(basis::SubBasis) = nvariables(parent(basis))
variables(basis::FullBasis) = variables(basis.map)
variables(basis::SubBasis) = variables(parent(basis))

function explicit_basis_covering(
    full::FullBasis{B},
    target::SubBasis{B},
) where {B}
    return SA.SubBasis(full, target.keys)
end

function monomial_type(::Type{<:SA.SubBasis{T,I,K,B}}) where {T,I,K,B}
    return monomial_type(B)
end

# FIXME workaround for TP, we should redirect to typeof
monomial_type(b::FullBasis) = monomial_type(b.map)
monomial_type(b::SubBasis) = monomial_type(parent(b))

# FIXME type piracy
SA.comparable(::ExponentsIterator{M}) where {M} = M()

function unsafe_basis(
    ::Type{B},
    monomials_vec::AbstractVector{<:AbstractMonomial},
) where {B<:AbstractMonomialIndexed}
    vars = variables(monomials_vec)
    basis = FullBasis{B}(vars)
    keys =
        is_commutative(vars) ?
        map(m -> _aligned_exponents(m, vars), monomials_vec) : monomials_vec
    return SA.SubBasis(basis, keys)
end

function Base.getindex(
    basis::FullBasis{B},
    monomials_vec::AbstractVector{M},
) where {B,M<:AbstractMonomial}
    return SA.SubBasis(
        basis,
        map(
            m -> exponents(m, variables(basis)),
            monomial_vector(monomials_vec),
        ),
    )
end

function SubBasis{B}(
    monomials_vec::AbstractVector{<:AbstractTermLike},
) where {B<:AbstractMonomialIndexed}
    return unsafe_basis(
        B,
        monomial_vector(monomials_vec)::AbstractVector{<:AbstractMonomial},
    )
end

SubBasis{B}(monos::Tuple) where {B} = SubBasis{B}([monos...])

function implicit_basis(basis::SA.SubBasis)
    return parent(basis)
end
function implicit_basis(basis::SA.ImplicitBasis)
    return basis
end

function implicit(a::SA.AlgebraElement)
    basis = implicit_basis(SA.basis(a))
    basis === SA.basis(a) && return a
    return algebra_element(SA.coeffs(a, basis), basis)
end

function implicit(t::SA.Term)
    basis = implicit_basis(SA.basis(parent(t)))
    basis === SA.basis(parent(t)) && return t
    return SA.Term(algebra(basis), basis[SA.basis_element(t)], coefficient(t))
end

function _comparable_type(
    ::Type{<:SA.MappedBasis{T,I,<:ExponentsIterator{M}}},
) where {T,I,M}
    return M
end
function _coeffs_type(
    ::Type{C},
    ::Type{B},
) where {C,B<:Union{FullBasis,_NoncommutativeBasis}}
    T = eltype(C) # Even works for `NTuple`!
    E = SA.key_type(B)
    L = _comparable_type(B)
    return SA.SparseCoefficients{E,T,_similar_type(C, E),C,L}
end

function _coeffs_type(::Type{C}, ::Type{B}) where {C,B<:SubBasis}
    return C
end

function algebra_element_type(::Type{C}, ::Type{B}) where {C,B}
    A = MA.promote_operation(algebra, B)
    return SA.AlgebraElement{eltype(C),A,_coeffs_type(C, B)}
end

function MA.promote_operation(
    ::typeof(implicit),
    ::Type{AE},
) where {AG,T,AE<:SA.AlgebraElement{T,AG}}
    BT =
        MA.promote_operation(implicit_basis, MA.promote_operation(SA.basis, AE))
    return algebra_element_type(Vector{T}, BT)
end

MA.promote_operation(::typeof(implicit_basis), B::Type{<:SA.ImplicitBasis}) = B
function MA.promote_operation(
    ::typeof(implicit_basis),
    ::Type{<:SA.SubBasis{T,I,K,B}},
) where {T,I,K,B}
    return B
end

function _explicit_basis(coeffs, basis::SA.ImplicitBasis)
    return SA.SubBasis(basis, _lazy_collect(SA.keys(coeffs)))
end

_explicit_basis(_, basis::SubBasis) = basis

function explicit_basis(p::AbstractTermLike)
    return explicit_basis(algebra_element(p))
end

function explicit_basis(a::SA.AlgebraElement)
    return _explicit_basis(SA.coeffs(a), SA.basis(a))
end

function explicit_basis_type(BT::Type{<:FullBasis{B,V,E}}) where {B,V,E}
    return SA.SubBasis{eltype(BT),Int,SA.key_type(BT),BT,Vector{E}}
end

function empty_basis(
    ::Type{<:SubBasis{B,M}},
) where {B<:AbstractMonomialIndexed,M}
    return unsafe_basis(B, empty_monomial_vector(M))
end

function maxdegree_basis(
    basis::FullBasis{B},
    maxdeg::Int,
) where {B<:AbstractMonomialIndexed}
    return unsafe_basis(B, monomials(variables(basis), 0:maxdeg))
end

variables(c::SA.AbstractCoefficients) = variables(SA.keys(c))

_lazy_collect(v::AbstractVector) = collect(v)
_lazy_collect(v::Tuple) = collect(v)
_lazy_collect(v::Vector) = v

sparse_coefficients(p::AbstractPolynomial) = SA.coeffs(p)
sparse_coefficients(t::AbstractTermLike) = SA.coeffs(algebra_element(t))
algebra_element(p::AbstractPolynomial) = p
algebra_element(t::AbstractTerm) = SA.algebra_element(t)
algebra_element(m::AbstractMonomialLike) = SA.algebra_element(term(m))

function algebra_element(f::Function, basis::SubBasis)
    return algebra_element(map(f, eachindex(basis)), basis)
end

_similar_type(::Type{<:NTuple{N,Any}}, ::Type{T}) where {N,T} = NTuple{N,T}
function _similar_type(::Type{V}, ::Type{T}) where {V<:AbstractVector,T}
    return SA.similar_type(V, T)
end

function full_basis_type(::Type{B}, ::Type{Polynomial{B,V,E}}) where {B,V,E}
    if !is_commutative(V)
        B === Monomial ||
            error("Noncommutative variables require the Monomial basis")
        return _NoncommutativeBasis{Polynomial{B,V,E},V}
    end
    return SA.MappedBasis{
        Polynomial{B,V,E},
        E,
        ExponentsIterator{ordering(V),Nothing,E},
        Variables{B,V},
        typeof(exponents),
    }
end

function MA.promote_operation(
    ::typeof(algebra_element),
    P::Type{<:AbstractPolynomialLike},
)
    return polynomial_type(P)
end

_one_if_type(α) = α
_one_if_type(::Type{T}) where {T} = one(T)

function constant_algebra_element_type(
    ::Type{BT},
    ::Type{T},
) where {BT<:Union{FullBasis,_NoncommutativeBasis},T}
    return algebra_element_type(Vector{T}, BT)
end

function constant_algebra_element(b::Union{FullBasis,_NoncommutativeBasis}, α)
    A = algebra(b)
    return SA.algebra_element(SA.Term(A, b[one(SA.object(A))], _one_if_type(α)))
end

function constant_algebra_element_type(
    ::Type{B},
    ::Type{T},
) where {B<:SubBasis,T}
    A = MA.promote_operation(algebra, B)
    return SA.AlgebraElement{T,A,Vector{T}}
end

function constant_algebra_element(basis::SubBasis, α)
    return algebra_element(
        [_one_if_type(α)],
        SA.SubBasis(
            parent(basis),
            [constant_monomial_exponents(parent(basis))],
        ),
    )
end

function SA.promote_with_map(::Variables{B}, vars, map) where {B}
    return Variables{B}(vars), map
end

SA.promote_with_map(::FullBasis, vars, m) = FullBasis(vars), m

function promote_variables_with_maps(
    a::Variables{B1},
    b::Variables{B2},
) where {B1,B2}
    B1 === B2 || throw(
        ArgumentError(
            "Cannot promote different basis families $B1 and $B2; use convert_basis explicitly.",
        ),
    )
    _a, _b = promote_variables_with_maps(a.variables, b.variables)
    return SA.maybe_promote(a, _a...), SA.maybe_promote(b, _b...)
end

function SA.promote_bases_with_maps(a::FullBasis, b::FullBasis)
    _a, _b = promote_variables_with_maps(a.map, b.map)
    return SA.maybe_promote(a, _a...), SA.maybe_promote(b, _b...)
end

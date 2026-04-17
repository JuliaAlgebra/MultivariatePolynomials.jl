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

struct _ExponentsWithVariables{V}
    variables::V
end

(e::_ExponentsWithVariables)(mono) = exponents(mono, e.variables)

function Base.:(==)(a::_ExponentsWithVariables, b::_ExponentsWithVariables)
    return unique(sort(a.variables)) == unique(sort(b.variables))
end

function FullBasis(vars::Variables{B,V}) where {B,V}
    O = ordering(vars.variables)
    exps = ExponentsIterator{O}(constant_monomial_exponents(vars))
    inverse_map = if is_commutative(V)
        exponents
    else
        _ExponentsWithVariables(vars.variables)
    end
    return SA.MappedBasis{Polynomial{B,V,eltype(exps)}}(exps, vars, inverse_map)
end

function FullBasis{B}(vars) where {B}
    return FullBasis(Variables{B,typeof(vars)}(vars))
end

function FullBasis{B}(p::_APL) where {B}
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
    return SA.SubBasis(
        FullBasis{B}(variables(monomials_vec)),
        exponents.(monomials_vec),
    )
end

function Base.getindex(
    ::FullBasis{B},
    monomials_vec::AbstractVector{M},
) where {B,M<:AbstractMonomial}
    return unsafe_basis(B, monomial_vector(monomials_vec)::AbstractVector{M})
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

implicit_basis(basis::SubBasis) = parent(basis)
implicit_basis(basis::FullBasis) = basis

function implicit(a::SA.AlgebraElement)
    basis = implicit_basis(SA.basis(a))
    return algebra_element(SA.coeffs(a, basis), basis)
end

function _coeffs_type(::Type{C}, ::Type{B}) where {C,B<:FullBasis}
    T = eltype(C)
    E = SA.key_type(B)
    return SA.SparseCoefficients{E,T,_similar_type(C, E),C,typeof(isless)}
end

function _coeffs_type(::Type{C}, ::Type{B}) where {C,B<:SubBasis}
    return C
end

function algebra_element_type(
    ::Type{C},
    ::Type{B},
) where {C,B}
    A = MA.promote_operation(algebra, B)
    return SA.AlgebraElement{
        eltype(C),
        A,
        _coeffs_type(C, B),
    }
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

function _explicit_basis(coeffs, basis::FullBasis{B}) where {B}
    return SA.SubBasis(basis, _lazy_collect(SA.keys(coeffs)))
end

_explicit_basis(_, basis::SubBasis) = basis

function explicit_basis(p::AbstractPolynomialLike)
    return SubBasis{Monomial}(monomials(p))
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

function sparse_coefficients(p::AbstractPolynomial)
    return SA.SparseCoefficients(
        exponents.(monomials(p)),
        _lazy_collect(coefficients(p)),
    )
end

function sparse_coefficients(t::AbstractTermLike)
    return SA.SparseCoefficients((exponents(t),), (coefficient(t),))
end

function algebra_element(p::_APL)
    return algebra_element(sparse_coefficients(p), FullBasis{Monomial}(p))
end

function algebra_element(f::Function, basis::SubBasis)
    return algebra_element(map(f, eachindex(basis)), basis)
end

_similar_type(::Type{<:NTuple{N,Any}}, ::Type{T}) where {N,T} = NTuple{N,T}
function _similar_type(::Type{V}, ::Type{T}) where {V<:AbstractVector,T}
    return SA.similar_type(V, T)
end

function full_basis_type(
    ::Type{B},
    ::Type{P},
) where {B,P<:AbstractPolynomialLike}
    V = MA.promote_operation(variables, P)
    E = _similar_type(V, Int)
    O = ordering(P)
    return SA.MappedBasis{
        Polynomial{B,V,E},
        E,
        ExponentsIterator{O,Nothing,E},
        Variables{B,V},
        typeof(exponents),
    }
end

function MA.promote_operation(
    ::typeof(algebra_element),
    P::Type{<:AbstractPolynomialLike{T}},
) where {T}
    return algebra_element_type(Vector{T}, full_basis_type(Monomial, P))
end

_one_if_type(α) = α
_one_if_type(::Type{T}) where {T} = one(T)

function constant_algebra_element_type(
    ::Type{BT},
    ::Type{T},
) where {BT<:FullBasis,T}
    return algebra_element_type(Tuple{T}, BT)
end

function constant_algebra_element(b::FullBasis, α)
    return algebra_element(
        SA.SparseCoefficients(
            (constant_monomial_exponents(b),),
            (_one_if_type(α),),
        ),
        b,
    )
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

function promote_variables_with_maps(a::Variables, b::Variables)
    _a, _b = promote_variables_with_maps(a.variables, b.variables)
    return SA.maybe_promote(a, _a...), SA.maybe_promote(b, _b...)
end

function SA.promote_bases_with_maps(a::FullBasis, b::FullBasis)
    _a, _b = promote_variables_with_maps(a.map, b.map)
    return SA.maybe_promote(a, _a...), SA.maybe_promote(b, _b...)
end

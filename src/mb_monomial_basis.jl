# From MultivariateBases/src/monomial.jl
# Monomial basis type and operations

# Note: AbstractMonomial is already defined in MP as:
#   abstract type AbstractMonomial <: AbstractMonomialLike end
# MB has its own AbstractMonomial <: AbstractMonomialIndexed.
# We use a different name to avoid collision.
abstract type AbstractMonomialBasis <: AbstractMonomialIndexed end

function explicit_basis_covering(
    full::FullBasis{B},
    target::SubBasis{<:AbstractMonomialBasis},
) where {B<:AbstractMonomialBasis}
    return SA.SubBasis(full, target.keys)
end

function explicit_basis_covering(
    full::FullBasis{B},
    target::SubBasis{B},
) where {B<:AbstractMonomialBasis}
    @assert full == parent(target)
    return SA.SubBasis(parent(target), target.keys)
end

function Base.adjoint(p::Polynomial{B}) where {B<:AbstractMonomialIndexed}
    mono = adjoint(monomial(p))
    return Polynomial(Variables{B}(variables(mono)), exponents(mono))
end

"""
    struct Monomial <: AbstractMonomialBasis end

Monomial basis with the monomials of the vector `monomials`.
For instance, `SubBasis{Monomial}([1, x, y, x^2, x*y, y^2])` is the monomial basis
for the subspace of quadratic polynomials in the variables `x`, `y`.
"""
struct Monomial <: AbstractMonomialBasis end

degree_one_univariate_polynomial(::Type{Monomial}, value) = value

function recurrence_eval(::Type{Monomial}, previous, value, degree)
    return previous[degree] * value
end

# /!\ assumes commutative variables
function (m::MStruct{Monomial,V,E})(a::E, b::E, ::Type{E}) where {V,E}
    return SA.SparseCoefficients((a .+ b,), (1,))
end

# Called by DiracMStructure for noncommutative variables
function Base.:*(a::Polynomial{Monomial}, b::Polynomial{Monomial})
    return Polynomial{Monomial}(monomial(a) * monomial(b))
end

SA.coeffs(p::Polynomial{Monomial}, ::FullBasis{Monomial}) = p.monomial

function polynomial_type(
    ::Union{SubBasis{B,V,E},Type{<:SubBasis{B,V,E}}},
    ::Type{T},
) where {B,V,E,T}
    return _polynomial_type(B, V, T)
end

function polynomial_type(
    ::Type{Polynomial{B,V,E}},
    ::Type{T},
) where {B,V,E,T}
    return _polynomial_type(B, V, T)
end

function keys_as_monomials(keys, mb::FullBasis)
    return monomial_vector(monomial.(Ref(mb.map.variables), keys))
end
keys_as_monomials(mb::SubBasis) = keys_as_monomials(mb.keys, parent(mb))

function polynomial(f::Function, mb::SubBasis{Monomial})
    return polynomial(f, keys_as_monomials(mb))
end

function polynomial(Q::AbstractMatrix, mb::SubBasis{Monomial}, T::Type)
    return polynomial(Q, keys_as_monomials(mb), T)
end

function coefficients(
    p::AbstractPolynomialLike,
    basis::SubBasis{Monomial},
)
    return coefficients(p, keys_as_monomials(basis))
end

function coefficients(p::AbstractPolynomialLike, ::FullBasis{Monomial})
    return p
end

function _assert_constant(α) end

function _assert_constant(
    x::Union{Polynomial,SA.AlgebraElement,AbstractPolynomialLike},
)
    return error("Expected constant element, got type `$(typeof(x))`")
end

mindegree(basis::SubBasis{Monomial}) = minimum(sum, basis.keys)
maxdegree(basis::SubBasis) = maximum(sum, basis.keys)
extdegree(basis::SubBasis{Monomial}) = extrema(sum, basis.keys)

variable_index(basis::FullBasis, v) = variable_index(basis.map, v)
variable_index(basis::SubBasis, v) = variable_index(basis.parent_basis, v)

function mindegree(basis::SubBasis{Monomial}, v)
    i = variable_index(basis, v)
    return minimum(Base.Fix2(getindex, i), basis.keys)
end
function maxdegree(basis::SubBasis, v)
    i = variable_index(basis, v)
    return maximum(Base.Fix2(getindex, i), basis.keys)
end
function extdegree(basis::SubBasis{Monomial}, v)
    i = variable_index(basis, v)
    return extrema(Base.Fix2(getindex, i), basis.keys)
end

_promote_coef(::Type{T}, ::Type{Monomial}) where {T} = T

function polynomial_type(::Type{<:FullBasis{B,V}}, ::Type{T}) where {T,B,V}
    return _polynomial_type(B, V, T)
end

function _polynomial_type(::Type{B}, ::Type{V}, ::Type{T}) where {B,V,T}
    return polynomial_type(monomial_type(V), _promote_coef(T, B))
end

_vec(v::Vector) = v
_vec(v::AbstractVector) = collect(v)

function SA.coeffs(
    cfs,
    source::MonomialIndexedBasis{B1},
    target::MonomialIndexedBasis{B2},
) where {B1,B2}
    source === target && return cfs
    source == target && return cfs
    if B1 === B2 && target isa FullBasis
        return SA.SparseCoefficients(_vec(source.keys), _vec(cfs))
    else
        res = SA.zero_coeffs(
            _promote_coef(_promote_coef(SA.value_type(cfs), B1), B2),
            target,
        )
        return SA.coeffs!(res, cfs, source, target)
    end
end

SA.star(::SubBasis, coeffs) = SA.star.(coeffs)

function _show_vector(io::IO, mime::MIME, v, map = identity)
    print(io, '[')
    first = true
    for el in v
        if !first
            print(io, ", ")
        end
        first = false
        show(io, mime, map(el))
    end
    return print(io, ']')
end

function _show(io::IO, mime::MIME, basis::SubBasis{B}) where {B}
    print(io, "SubBasis{$(nameof(B))}(")
    _show_vector(
        io,
        mime,
        basis.keys,
        Base.Fix1(monomial, variables(basis)),
    )
    print(io, ')')
    return
end

function _show(io::IO, mime::MIME, basis::FullBasis{B}) where {B}
    print(io, "FullBasis{$(nameof(B))}(")
    _show_vector(io, mime, variables(basis))
    print(io, ")")
    return
end

function Base.show(
    io::IO,
    mime::MIME"text/plain",
    basis::Union{Variables,SubBasis,FullBasis},
)
    return _show(io, mime, basis)
end

function Base.show(
    io::IO,
    mime::MIME"text/print",
    basis::Union{Variables,SubBasis,FullBasis},
)
    return _show(io, mime, basis)
end

function Base.print(io::IO, basis::Union{Variables,SubBasis,FullBasis})
    return show(io, MIME"text/print"(), basis)
end

function Base.show(io::IO, basis::Union{Variables,SubBasis,FullBasis})
    return show(io, MIME"text/plain"(), basis)
end

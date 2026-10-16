# Polynomial{Monomial,...} acts as a monomial: support term creation and const mult
term(coeff, p::Polynomial{<:AbstractMonomialBasis}) = SA.Term(coeff, p)
term(p::Polynomial{<:AbstractMonomialBasis}) = SA.Term(one(Int), p)
left_constant_mult(α, p::Polynomial{<:AbstractMonomialBasis}) = SA.Term(α, p)
function right_constant_mult(p::Polynomial{<:AbstractMonomialBasis}, α)
    return left_constant_mult(α, p)
end

function term_type(::Type{Polynomial{B,V,E}}, ::Type{T}) where {B,V,E,T}
    BT = full_basis_type(B, Polynomial{B,V,E})
    A = MA.promote_operation(algebra, BT)
    return SA.Term{T,A,SA.key_type(BT)}
end
term_type(::Type{P}) where {P<:Polynomial} = term_type(P, Int)

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

degree_one_univariate_polynomial(::Type{Monomial}, value) = value

function recurrence_eval(::Type{Monomial}, previous, value, degree)
    return previous[degree] * value
end

# /!\ assumes commutative variables
function (m::MStruct{Monomial,V,E})(a::E, b::E, ::Type{E}) where {V,E}
    return SA.SparseCoefficients((a .+ b,), (1,))
end

function SA.term_product_style(
    ::MStruct{Monomial,V,E,E,<:FullBasis},
    ::AbstractMonomialOrdering,
) where {V,E}
    return SA.OrderedTermProduct()
end

# Monomial multiplication: align variables via promote_variables_with_maps,
# then apply f to the aligned exponent vectors.
function map_exponents(f, a::Polynomial{Monomial}, b::Polynomial{Monomial})
    a, b = SA.promote_bases(a, b)
    return Polynomial(a.variables, f.(exponents(a), exponents(b)))
end

function _aligned_exponents(m::AbstractMonomialLike, vars)
    m = monomial(m)
    return variables(m) == vars ? exponents(m) : exponents(m, vars)
end

function map_exponents_to!(
    output::Polynomial{Monomial,V,<:AbstractVector},
    f::F,
    a::AbstractMonomialLike,
    b::AbstractMonomialLike,
) where {V,F}
    vars = variables(output)
    ea = _aligned_exponents(a, vars)
    eb = _aligned_exponents(b, vars)
    broadcast!(f, exponents(output), ea, eb)
    return output
end

function map_exponents_to!(
    output::Polynomial{Monomial,V,<:Tuple},
    f::F,
    a::AbstractMonomialLike,
    b::AbstractMonomialLike,
) where {V,F}
    vars = variables(output)
    ea = _aligned_exponents(a, vars)
    eb = _aligned_exponents(b, vars)
    return Polynomial(output.variables, map(f, ea, eb))
end

function map_exponents!(f, a::AbstractMonomial, b::AbstractMonomialLike)
    return map_exponents_to!(a, f, a, b)
end

function Base.:*(a::Polynomial{Monomial}, b::Polynomial{Monomial})
    return map_exponents(+, a, b)
end

# Monomial power: multiply exponents by n
function Base.:^(m::Polynomial{Monomial}, n::Integer)
    return Polynomial(m.variables, m.exponents .* n)
end

SA.coeffs(p::Polynomial{Monomial}, ::FullBasis{Monomial}) = p

function polynomial_type(
    ::Union{SubBasis{B,V,E},Type{<:SubBasis{B,V,E}}},
    ::Type{T},
) where {B,V,E,T}
    return polynomial_type(Polynomial{B,V,E}, T)
end

function polynomial_type(::Type{P}, ::Type{T}) where {B,P<:Polynomial{B},T}
    return algebra_element_type(Vector{T}, full_basis_type(B, P))
end

function keys_as_monomials(keys, mb::FullBasis)
    return monomial_vector(monomial.(Ref(mb.map.variables), keys))
end
keys_as_monomials(mb::SubBasis) = keys_as_monomials(mb.keys, parent(mb))

function polynomial(f::Function, mb::SubBasis)
    return implicit(algebra_element(f, mb))
end

function polynomial(Q::AbstractMatrix, mb::SubBasis{Monomial}, T::Type)
    return polynomial(Q, keys_as_monomials(mb), T)
end

function coefficients(p::AbstractPolynomialLike, basis::SubBasis{Monomial})
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

mindegree(basis::SubBasis{Monomial}) = minimum(sum ∘ exponents, basis)
maxdegree(basis::SubBasis) = maximum(sum ∘ exponents, basis)
extdegree(basis::SubBasis{Monomial}) = extrema(sum ∘ exponents, basis)

variable_index(basis::FullBasis, v) = variable_index(basis.map, v)
variable_index(basis::SubBasis, v) = variable_index(basis.parent_basis, v)

function mindegree(basis::SubBasis{Monomial}, v)
    return minimum(Base.Fix2(degree, v) ∘ monomial, basis)
end
function maxdegree(basis::SubBasis, v)
    return maximum(Base.Fix2(degree, v) ∘ monomial, basis)
end
function extdegree(basis::SubBasis{Monomial}, v)
    return extrema(Base.Fix2(degree, v) ∘ monomial, basis)
end

_promote_coef(::Type{T}, ::Type{Monomial}) where {T} = T

function polynomial_type(::Type{BT}, ::Type{T}) where {T,BT<:FullBasis}
    return algebra_element_type(Vector{T}, BT)
end

function _polynomial_type(::Type{B}, ::Type{V}, ::Type{T}) where {B,V,T}
    return polynomial_type(monomial_type(V), _promote_coef(T, B))
end

_vec(v::Vector) = v
_vec(v::AbstractVector) = collect(v)

function SA.coeffs(cfs, source::SubBasis{B}, target::FullBasis{B}) where {B}
    return SA.SparseCoefficients(
        _vec(source.keys),
        _vec(cfs),
        SA.comparable(target),
    )
end

function SA.promote_coefficient_type(
    ::Union{typeof(SA.coeffs),typeof(SA.adjoint_coeffs)},
    ::Type{T},
    ::Type{<:MonomialIndexedBasis{B1}},
    ::Type{<:MonomialIndexedBasis{B2}},
) where {T,B1,B2}
    return B1 === B2 ? T : _promote_coef(_promote_coef(T, B1), B2)
end

function SA.adjoint_coeffs!(
    res,
    cfs,
    source::MonomialIndexedBasis{B1},
    target::MonomialIndexedBasis{B2},
) where {B1,B2}
    if B1 === B2
        return @invoke SA.adjoint_coeffs!(
            res,
            cfs,
            source::SA.AbstractBasis,
            target::SA.AbstractBasis,
        )
    end
    MA.operate!(zero, res)
    target_full = implicit_basis(target)
    for (i, src_elem) in enumerate(source)
        full_col = SA.coeffs(algebra_element(src_elem), target_full)
        col = SA.coeffs(full_col, target_full, target)
        for (j, v) in SA.nonzero_pairs(col)
            res[i] += v * cfs[j]
        end
    end
    return res
end

# FIXME this assumes that the basis is invariant under adjoint
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
    _show_vector(io, mime, basis.keys, Base.Fix1(monomial, variables(basis)))
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

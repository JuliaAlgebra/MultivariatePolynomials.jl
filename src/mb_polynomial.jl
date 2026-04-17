# From MultivariateBases/src/polynomial.jl
# Polynomial{B,V,E} is a basis element (single monomial-like indexed object)

# Bridge code for SA types and MP types
variables(p::SA.AlgebraElement) = variables(explicit_basis(p))
Base.keytype(p::AbstractPolynomialLike) = monomial_type(p)
SA.value_type(p::AbstractPolynomialLike) = coefficient_type(p)
SA.nonzero_pairs(p::AbstractPolynomialLike) = terms(p)
function Base.similar(p::PT, ::Type{T}) where {PT<:AbstractPolynomial,T}
    return convert(similar_type(PT, T), copy(p))
end
function Base.getindex(p::AbstractPolynomialLike, mono::AbstractMonomial)
    return coefficient(p, mono)
end
Base.iterate(t::SA.Term) = iterate(t, 1)
function Base.iterate(t::SA.Term, state)
    if state == 1
        return monomial(t), 2
    elseif state == 2
        return coefficient(t), 3
    else
        return nothing
    end
end
function SA.unsafe_push!(p::AbstractPolynomial, mono::AbstractMonomial, α)
    return MA.operate!(MA.add_mul, p, α, mono)
end
function MA.operate!(
    ::SA.UnsafeAddMul{typeof(*)},
    mc::AbstractPolynomial,
    val,
    c::AbstractPolynomialLike,
)
    return MA.operate!(MA.add_mul, mc, val, c)
end
MA.operate!(::typeof(SA.canonical), p::AbstractPolynomial) = p
function MA.promote_operation(
    ::typeof(SA.canonical),
    ::Type{P},
) where {P<:AbstractPolynomialLike}
    return P
end

abstract type AbstractMonomialIndexed end

"""
    struct Polynomial{B<:AbstractMonomialIndexed,V,E}
        variables::Variables{B,V}
        exponents::E
    end

Polynomial of basis `FullBasis{B,V,E}(variables)` at index `exponents`.
This represents a single basis element, not a sum of terms.
"""
struct Polynomial{B<:AbstractMonomialIndexed,V,E}
    variables::Variables{B,V}
    exponents::E
end

function Polynomial(v::Variables{B,V}, e) where {B,V}
    return Polynomial{B,V,typeof(e)}(v, e)
end

function Polynomial{B}(v::AbstractVariable) where {B}
    vars = variables(v)
    return Polynomial(Variables{B}(vars), map(_ -> 1, vars))
end

function Polynomial{B}(mono::AbstractMonomial) where {B}
    vars = variables(mono)
    return Polynomial(Variables{B}(vars), exponents(mono))
end

exponents(p::Polynomial) = p.exponents
exponents(p::Polynomial, vars) = exponents(monomial(p), vars)
monomial(p::Polynomial) = monomial(variables(p), exponents(p))

function Base.hash(p::Polynomial{B}, u::UInt) where {B}
    return hash(p.variables, hash(p.exponents, u))
end

function Base.isequal(p::Polynomial{B}, q::Polynomial{B}) where {B}
    return isequal(p.variables, q.variables) &&
           isequal(p.exponents, q.exponents)
end

Base.isone(p::Polynomial) = all(iszero, p.exponents)
isconstant(p::Polynomial) = all(iszero, p.exponents)
# A monomial basis element is never zero
Base.iszero(p::Polynomial) = false
# constant_monomial: return a monomial with all-zero exponents
constant_monomial(p::Polynomial) = Polynomial(p.variables, zero(p.exponents))

Base.iterate(p::Polynomial) = p, nothing
Base.iterate(::Polynomial, ::Nothing) = nothing

function Base.:(==)(p::Polynomial{B}, q::Polynomial{B}) where {B}
    return p.variables == q.variables && p.exponents == q.exponents
end

variables(p::Polynomial) = variables(p.variables)
nvariables(p::Polynomial) = nvariables(p.variables)

monomial_type(::Type{<:SA.SparseCoefficients{K}}) where {K} = K
# SA.Term constructor for Polynomial{...} basis elements is in mb_monomial_basis.jl

polynomial(p::Polynomial) = polynomial(algebra_element(p))

function algebra_element(p, basis::SA.AbstractBasis)
    return SA.AlgebraElement(p, algebra(basis))
end

function _algebra_element(p, ::Type{B}) where {B<:AbstractMonomialIndexed}
    return algebra_element(
        sparse_coefficients(p),
        FullBasis{B}(variables(p)),
    )
end

function algebra_element(p::Polynomial{B}) where {B}
    return _algebra_element(monomial(p), B)
end

function Base.:*(a::Polynomial{B}, b::SA.AlgebraElement) where {B}
    return _algebra_element(a) * b
end

function _show(io::IO, mime::MIME, p::Polynomial{B}) where {B}
    if B != Monomial
        print(io, B)
        print(io, "(")
    end
    # Display monomial from variables and exponents directly
    vars = variables(p)
    exps = exponents(p)
    if all(iszero, exps)
        print(io, "1")
    else
        first = true
        for (v, e) in zip(vars, exps)
            iszero(e) && continue
            if !first
                print(io, "*")
            end
            first = false
            show(io, mime, v)
            if e > 1
                print(io, "^", e)
            end
        end
    end
    if B != Monomial
        print(io, ")")
    end
    return
end

function Base.show(io::IO, mime::MIME"text/latex", p::Polynomial)
    print(io, "\$\$ ")
    _show(io, mime, p)
    print(io, " \$\$")
    return
end

function Base.show(io::IO, mime::MIME"text/plain", p::Polynomial)
    return _show(io, mime, p)
end

function Base.show(io::IO, mime::MIME"text/print", p::Polynomial)
    return _show(io, mime, p)
end

Base.show(io::IO, p::Polynomial) = show(io, MIME"text/plain"(), p)
Base.print(io::IO, p::Polynomial) = show(io, MIME"text/print"(), p)

# zero for a Polynomial type: return an AlgebraElement with no terms
function Base.zero(::Type{Polynomial{B,V,E}}) where {B,V,E}
    vars = V()
    basis = FullBasis{B}(vars)
    sc = SA.SparseCoefficients(E[], Int[])
    return SA.algebra_element(sc, algebra(basis))
end

function Base.zero(p::Polynomial)
    basis = FullBasis{typeof_basis(p)}(p)
    sc = SA.SparseCoefficients(typeof(p.exponents)[], Int[])
    return SA.algebra_element(sc, algebra(basis))
end
typeof_basis(::Polynomial{B}) where {B} = B

function convert_basis(basis::SA.AbstractBasis, p::AbstractPolynomialLike)
    return convert_basis(basis, _algebra_element(p, Monomial))
end

function convert_basis(basis::SA.AbstractBasis, p::SA.AlgebraElement)
    return SA.AlgebraElement(SA.coeffs(p, basis), algebra(basis))
end

# isapprox: AlgebraElement IS the polynomial, so compare via algebra_element
function Base.isapprox(
    p::AbstractPolynomialLike,
    a::SA.AlgebraElement;
    kws...,
)
    return isapprox(algebra_element(p), a; kws...)
end

function Base.isapprox(
    a::SA.AlgebraElement,
    p::AbstractPolynomialLike;
    kws...,
)
    return isapprox(p, a; kws...)
end

function Base.isapprox(a::SA.AlgebraElement, α::Number; kws...)
    return isapprox(
        a,
        α * constant_algebra_element(SA.basis(a), typeof(α));
        kws...,
    )
end

# Type operations for AlgebraElement
function monomial_type(::Type{<:SA.AlgebraElement{T,A}}) where {T,A}
    return monomial_type(A)
end
function polynomial_type(::Type{<:SA.AlgebraElement{T,A}}) where {A,T}
    return polynomial_type(A, T)
end
monomial_type(::Type{<:SA.StarAlgebra{O}}) where {O} = monomial_type(O)
function polynomial_type(::Type{A}, ::Type{T}) where {A<:SA.StarAlgebra,T}
    return polynomial_type(MA.promote_operation(SA.basis, A), T)
end

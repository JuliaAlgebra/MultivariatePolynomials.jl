variables(p::AbstractPolynomial) = variables(SA.basis(p))

function BasisElement{B}(v::AbstractVariable) where {B}
    vars = variables(v)
    return BasisElement(Variables{B}(vars), map(_ -> 1, vars))
end

function BasisElement{B}(mono::AbstractMonomial) where {B}
    vars = variables(mono)
    return BasisElement(Variables{B}(vars), exponents(mono))
end

exponents(p::BasisElement) = p.exponents
function exponents(p::BasisElement, vars)
    (all_vars, map), _ = promote_variables_with_maps(variables(p), vars)
    if all_vars != vars
        throw(
            ArgumentError(
                "The supplied variables cannot represent this monomial",
            ),
        )
    end
    return map === nothing ? MA.copy_if_mutable(exponents(p)) :
           map(exponents(p))
end

function Base.copy(p::BasisElement{B}) where {B}
    return BasisElement(Variables{B}(copy(variables(p))), copy(exponents(p)))
end

function Base.isequal(p::BasisElement{B}, q::BasisElement{B}) where {B}
    @assert is_commutative(variables(p)) && is_commutative(variables(q))
    a, b = SA.promote_bases(p, q)
    return isequal(exponents(a), exponents(b))
end

Base.isone(p::BasisElement) = all(iszero, p.exponents)
# A monomial basis element is never zero
Base.iszero(p::BasisElement) = false
# Ordering: derive from the variable type
ordering(::Type{BasisElement{B,V,E}}) where {B,V,E} = ordering(V)
# constant_monomial: return a monomial with all-zero exponents
function constant_monomial(p::BasisElement)
    return BasisElement(p.variables, map(zero, p.exponents))
end

Base.iterate(p::BasisElement) = p, nothing
Base.iterate(::BasisElement, ::Nothing) = nothing

function Base.:(==)(p::BasisElement{B}, q::BasisElement{B}) where {B}
    @assert is_commutative(variables(p)) && is_commutative(variables(q))
    a, b = SA.promote_bases(p, q)
    return exponents(a) == exponents(b)
end

variables(p::BasisElement) = variables(p.variables)

monomial_type(::Type{<:SA.SparseCoefficients{K}}) where {K} = K
monomial_type(::Type{<:SA.StarAlgebra{O,P}}) where {O<:Variables,P} = P

"""
    convert_basis(basis, p)

Express `p` in `basis`. Polynomial construction preserves its basis family;
arithmetic between different families requires an explicit conversion first.
"""
function convert_basis(
    basis::SA.AbstractBasis,
    p::Union{AbstractTermLike,BasisElement},
)
    return convert_basis(basis, polynomial(p))
end

function convert_basis(basis::SA.AbstractBasis, p::SA.AlgebraElement)
    return SA.AlgebraElement(SA.coeffs(p, basis), algebra(basis))
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
monomial_type(::Type{<:SA.StarAlgebra{O}}) where {O} = monomial_type(O)
function polynomial_type(::Type{A}, ::Type{T}) where {A<:SA.StarAlgebra,T}
    return polynomial_type(MA.promote_operation(SA.basis, A), T)
end

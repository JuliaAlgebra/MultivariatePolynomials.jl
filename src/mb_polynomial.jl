variables(p::AbstractPolynomial) = variables(SA.basis(p))

function Polynomial{B}(v::AbstractVariable) where {B}
    vars = variables(v)
    return Polynomial(Variables{B}(vars), map(_ -> 1, vars))
end

function Polynomial{B}(mono::AbstractMonomial) where {B}
    vars = variables(mono)
    return Polynomial(Variables{B}(vars), exponents(mono))
end

exponents(p::Polynomial) = p.exponents
function exponents(p::Polynomial, vars)
    (all_vars, map), _ = promote_variables_with_maps(variables(p), vars)
    if all_vars != vars
        throw(
            ArgumentError("The supplied variables cannot represent this monomial"),
        )
    end
    return map === nothing ? copy(exponents(p)) : map(exponents(p))
end
monomial(p::Polynomial) = monomial(variables(p), exponents(p))
monomial(p::AbstractMonomial) = p

function Base.copy(p::Polynomial{B}) where {B}
    return Polynomial(Variables{B}(copy(variables(p))), copy(exponents(p)))
end

function Base.hash(p::Polynomial{B}, u::UInt) where {B}
    return hash(
        Tuple(
            (v, e) for (v, e) in zip(variables(p), exponents(p)) if !iszero(e)
        ),
        hash(B, u),
    )
end
function Base.isequal(p::Polynomial{B}, q::Polynomial{B}) where {B}
    a, b = SA.promote_bases(p, q)
    return isequal(exponents(a), exponents(b))
end

Base.isone(p::Polynomial) = all(iszero, p.exponents)
isconstant(p::Polynomial) = all(iszero, p.exponents)
# A monomial basis element is never zero
Base.iszero(p::Polynomial) = false
# Polynomial{Monomial,...} is its own monomial type
monomial_type(::Type{PT}) where {PT<:Polynomial{<:AbstractMonomialIndexed}} = PT
monomial_type(p::Polynomial{<:AbstractMonomialIndexed}) = typeof(p)
# Ordering: derive from the variable type
ordering(::Type{Polynomial{B,V,E}}) where {B,V,E} = ordering(V)
# constant_monomial: return a monomial with all-zero exponents
constant_monomial(p::Polynomial) = Polynomial(p.variables, zero(p.exponents))

Base.iterate(p::Polynomial) = p, nothing
Base.iterate(::Polynomial, ::Nothing) = nothing

function Base.:(==)(p::Polynomial{B}, q::Polynomial{B}) where {B}
    a, b = SA.promote_bases(p, q)
    return exponents(a) == exponents(b)
end

variables(p::Polynomial) = variables(p.variables)
nvariables(p::Polynomial) = nvariables(p.variables)

monomial_type(::Type{<:SA.SparseCoefficients{K}}) where {K} = K
monomial_type(::Type{<:SA.StarAlgebra{O,P}}) where {O<:Variables,P} = P
# SA.Term constructor for Polynomial{...} basis elements is in mb_monomial_basis.jl


function algebra_element(p, basis::SA.AbstractBasis)
    return SA.AlgebraElement(p, algebra(basis))
end

function _algebra_element(p, ::Type{B}) where {B<:AbstractMonomialIndexed}
    return algebra_element(
        sparse_coefficients(p),
        FullBasis{B}(variables(p)),
    )
end

algebra_element(p::Polynomial{B}) where {B} = _algebra_element(monomial(p), B)
polynomial(p::Polynomial) = algebra_element(p)

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

function convert_basis(basis::SA.AbstractBasis, p::AbstractTermLike)
    return convert_basis(basis, _algebra_element(p, Monomial))
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

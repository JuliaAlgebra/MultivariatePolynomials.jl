# Variables and monomials enter the shared term/algebra arithmetic.
for op in [:+, :-]
    @eval begin
        Base.$op(a::AbstractMonomialLike, b::AbstractTerm) = $op(term(a), b)
        Base.$op(a::AbstractTerm, b::AbstractMonomialLike) = $op(a, term(b))
        Base.$op(a::AbstractMonomialLike, b::AbstractMonomialLike) =
            $op(term(a), term(b))
    end
end
function Base.:*(a::AbstractMonomialLike, b::AbstractTerm)
    return term(MA.copy_if_mutable(coefficient(b)), a * monomial(b))
end
function Base.:*(a::AbstractTerm, b::AbstractMonomialLike)
    return term(MA.copy_if_mutable(coefficient(a)), monomial(a) * b)
end
function Base.:*(a::AbstractMonomialLike, b::AbstractMonomialLike)
    return monomial(a) * monomial(b)
end
for op in (:+, :-, :*)
    @eval begin
        Base.$op(p::AbstractPolynomialLike, s::LinearAlgebra.UniformScaling) =
            $op(p, s.λ)
        Base.$op(s::LinearAlgebra.UniformScaling, p::AbstractPolynomialLike) =
            $op(s.λ, p)
    end
end

_term(α, mono) = term(α, MA.copy_if_mutable(mono))

function Base.isapprox(t1::AbstractTermLike, t2::AbstractTermLike; kwargs...)
    return isapprox(coefficient(t1), coefficient(t2); kwargs...) &&
           monomial(t1) == monomial(t2)
end
# `MA.operate(-, p)` redirects to `-p` as it assumes that `-p` can be modified
# through the MA API without modifying `p`. We should either copy the monomial
# here or implement a `MA.operate(-, p)` that copies it. We choose the first
# option.
Base.:-(m::AbstractMonomialLike) = _term(-1, MA.copy_if_mutable(m))
Base.:-(t::AbstractTermLike) = _term(MA.operate(-, coefficient(t)), monomial(t))
Base.:+(p::Union{AbstractMonomialLike,RationalPoly}) = p
Base.:*(p::Union{AbstractMonomialLike,RationalPoly}) = p

# Coefficients and variables commute
left_constant_mult(α, m::_BasisElementLike) = term(α, m)
right_constant_mult(m::_BasisElementLike, α) = left_constant_mult(α, m)

function left_constant_mult(α, t::SA.Term)
    return term(α * coefficient(t), monomial(t))
end
function right_constant_mult(t::SA.Term, α)
    return term(coefficient(t) * α, monomial(t))
end

function left_constant_mult(α, p::AbstractPolynomial)
    return map_coefficients(Base.Fix1(*, α), p)
end
function right_constant_mult(p::AbstractPolynomial, α)
    return map_coefficients(Base.Fix2(*, α), p)
end

LinearAlgebra.adjoint(v::AbstractVariable) = conj(v)
function LinearAlgebra.adjoint(p::BasisElement{B}) where {B}
    mono = adjoint(monomial(p))
    return BasisElement(Variables{B}(variables(mono)), exponents(mono))
end
LinearAlgebra.adjoint(m::AbstractMonomial) = conj(m)
function LinearAlgebra.adjoint(t::AbstractTerm)
    return _term(adjoint(coefficient(t)), adjoint(monomial(t)))
end
function LinearAlgebra.adjoint(p::AbstractPolynomial)
    return polynomial(map(LinearAlgebra.adjoint, terms(p)))
end
function LinearAlgebra.adjoint(r::RationalPoly)
    return adjoint(numerator(r)) / adjoint(denominator(r))
end
LinearAlgebra.hermitian_type(::Type{T}) where {T<:AbstractPolynomialLike} = T
function LinearAlgebra.hermitian(v::AbstractPolynomialLike, ::Symbol)
    isreal(v) || error(
        "Complex-valued polynomials cannot be interpreted as hermitian scalars",
    )
    return v
end
# This is the same implementation as `LinearAlgebra.ishermitian(::Number)`
LinearAlgebra.ishermitian(p::AbstractPolynomialLike) = p == conj(p)

LinearAlgebra.transpose(v::AbstractVariable) = v
LinearAlgebra.transpose(m::AbstractMonomial) = m
function LinearAlgebra.transpose(t::AbstractTerm)
    return _term(LinearAlgebra.transpose(coefficient(t)), monomial(t))
end
function LinearAlgebra.transpose(p::AbstractPolynomial)
    return map_coefficients(LinearAlgebra.transpose, p; nonzero = true)
end
function LinearAlgebra.transpose(r::RationalPoly)
    return transpose(numerator(r)) / transpose(denominator(r))
end

function LinearAlgebra.dot(p1::AbstractMonomialLike, p2::AbstractPolynomialLike)
    return p1' * p2
end
function LinearAlgebra.dot(
    p1::Union{AbstractPolynomial,AbstractTerm},
    p2::AbstractMonomialLike,
)
    return p1' * p2
end

_sum_product_operand(p) = p
_sum_product_operand(m::AbstractMonomialLike) = term(m)

function MA.operate!(
    op::Union{typeof(+),typeof(-)},
    p::AbstractPolynomial,
    s::LinearAlgebra.UniformScaling,
)
    return MA.operate!(op, p, s.λ)
end

function MA.operate!(
    op::Union{typeof(+),typeof(-)},
    p::AbstractPolynomial,
    m::AbstractMonomialLike,
)
    return MA.operate!(op, p, term(m))
end

function MA.operate_to!(
    output::AbstractPolynomial,
    op::Union{typeof(+),typeof(-),typeof(*)},
    p::AbstractPolynomial,
    m::AbstractMonomialLike,
)
    return MA.operate_to!(output, op, p, term(m))
end

function MA.operate_to!(
    output::AbstractPolynomial,
    ::typeof(*),
    m::AbstractMonomialLike,
    p::AbstractPolynomial,
)
    return MA.operate_to!(output, *, term(m), p)
end

for (A, B) in (
    (AbstractMonomialLike, AbstractPolynomialLike),
    (Union{AbstractPolynomial,AbstractTerm}, AbstractMonomialLike),
)
    @eval function MA.operate!(
        op::MA.AddSubMul,
        p::AbstractPolynomial,
        a::$A,
        b::$B,
    )
        return MA.operate!(
            op,
            p,
            _sum_product_operand(a),
            _sum_product_operand(b),
        )
    end
end

for (A, B) in (
    (AbstractMonomialLike, AbstractPolynomialLike),
    (Union{AbstractPolynomial,AbstractTerm}, AbstractMonomialLike),
)
    @eval function MA.operate(
        ::typeof(LinearAlgebra.dot),
        a::AbstractArray{<:$A},
        b::AbstractArray{<:$B},
    )
        # Conjugation can change the basis, so do it before basis promotion.
        return LinearAlgebra._dot_nonrecursive(adjoint.(a), b)
    end
end

for (A, B) in (
    (AbstractMonomialLike, AbstractPolynomialLike),
    (Union{AbstractPolynomial,AbstractTerm}, AbstractMonomialLike),
)
    @eval function LinearAlgebra._dot_nonrecursive(
        a::AbstractArray{<:$A},
        b::AbstractArray{<:$B},
    )
        return MA.fused_map_reduce(
            MA.add_mul,
            map(_sum_product_operand, a),
            map(_sum_product_operand, b),
        )
    end
end

# A bare monomial takes its coefficient directly from the scalar factor.
# Multiplying an existing term also multiplies its coefficient.
function _sum_monomial_products(a, b)
    MA._check_same_length(a, b)
    return sum(map(term, vec(a), vec(b)))
end

LinearAlgebra.symmetric_type(PT::Type{<:_APL}) = PT
LinearAlgebra.symmetric(p::_APL, ::Symbol) = p
LinearAlgebra.issymmetric(::_APL) = true

# Amazingly, this works! Thanks, StaticArrays.jl!
"""
Convert a tuple of variables into a static vector to allow array-like usage.
The element type of the vector will be Monomial{vars, length(vars)}.
"""
Base.vec(vars::Tuple{Vararg{AbstractVariable}}) = [vars...]
# vec(vars::Tuple{Vararg{TypedVariable}}) = SVector(vars)

# https://github.com/JuliaLang/julia/pull/23332
Base.:^(x::AbstractMonomialLike, p::Integer) = Base.power_by_squaring(x, p)
# Powers of SA.Term and SA.AlgebraElement are defined in StarAlgebras.

function MA.operate_to!(output, ::typeof(left_constant_mult), α, p::_APL)
    return SA.map_coefficients_to!(output, Base.Fix1(*, α), p)
end
function MA.operate_to!(output, ::typeof(right_constant_mult), p::_APL, α)
    return SA.map_coefficients_to!(output, Base.Fix2(*, α), p)
end
function MA.operate!(::typeof(left_constant_mult), α, p::_APL)
    return SA.map_coefficients!(Base.Fix1(*, α), p)
end
function MA.operate!(::typeof(right_constant_mult), p::_APL, α)
    return SA.map_coefficients!(Base.Fix2(MA.mul!!, α), p)
end

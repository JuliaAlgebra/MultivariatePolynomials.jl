# Variables and monomials enter the shared term/algebra arithmetic.
for op in [:+, :-, :*]
    @eval begin
        Base.$op(a::AbstractMonomialLike, b::AbstractTerm) = $op(term(a), b)
        Base.$op(a::AbstractTerm, b::AbstractMonomialLike) = $op(a, term(b))
    end
end
for op in [:+, :-]
    @eval Base.$op(a::AbstractMonomialLike, b::AbstractMonomialLike) =
        $op(term(a), term(b))
end
function Base.:*(a::AbstractMonomialLike, b::AbstractMonomialLike)
    return monomial(a) * monomial(b)
end
Base.:*(a::Number, b::AbstractMonomialLike) = term(a, monomial(b))
Base.:*(a::AbstractMonomialLike, b::Number) = term(b, monomial(a))

for op in [:+, :-]
    @eval begin
        Base.$op(p::AbstractMonomialLike, a::Number) =
            $op(algebra_element(p), a)
        Base.$op(a::Number, p::AbstractMonomialLike) =
            $op(a, algebra_element(p))
        Base.$op(p::AbstractTerm{T}, a::Union{T,Number}) where {T} =
            $op(algebra_element(p), a)
        Base.$op(a::Union{T,Number}, p::AbstractTerm{T}) where {T} =
            $op(a, algebra_element(p))
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
Base.:+(p::Union{_APL,RationalPoly}) = p
Base.:*(p::Union{_APL,RationalPoly}) = p

# Coefficients and variables commute
left_constant_mult(α, v::AbstractMonomialLike) = SA.Term(α, monomial(v))
right_constant_mult(m::AbstractMonomialLike, α) = left_constant_mult(α, m)
# Polynomial{Monomial,...} methods added in mb_monomial_basis.jl

function left_constant_mult(α, t::SA.Term)
    return term(α * coefficient(t), monomial(t))
end
function right_constant_mult(t::SA.Term, α)
    return term(coefficient(t) * α, monomial(t))
end

function left_constant_mult(α, p::AbstractPolynomialLike)
    return map_coefficients(Base.Fix1(*, α), p)
end
function right_constant_mult(p::AbstractPolynomialLike, α)
    return map_coefficients(Base.Fix2(*, α), p)
end

LinearAlgebra.adjoint(v::AbstractVariable) = conj(v)
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

function LinearAlgebra.dot(
    p1::AbstractPolynomialLike,
    p2::AbstractPolynomialLike,
)
    return p1' * p2
end
LinearAlgebra.dot(x::Number, p::AbstractPolynomialLike) = x' * p
LinearAlgebra.dot(p::AbstractPolynomialLike, x::Number) = p' * x

_sum_product_operand(p::Union{AbstractPolynomial,AbstractTerm}) = p
_sum_product_operand(m::AbstractMonomialLike) = term(m)

for (A, B) in (
    (AbstractPolynomialLike, AbstractPolynomialLike),
    (Number, AbstractPolynomialLike),
    (AbstractPolynomialLike, Number),
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

function LinearAlgebra._dot_nonrecursive(
    a::AbstractArray{<:AbstractPolynomialLike},
    b::AbstractArray{<:AbstractPolynomialLike},
)
    return SA.sum_products(
        map(_sum_product_operand, a),
        map(_sum_product_operand, b),
    )
end

function LinearAlgebra._dot_nonrecursive(
    a::AbstractArray{<:Number},
    b::AbstractArray{<:Union{AbstractPolynomial,AbstractTerm}},
)
    return SA.sum_products(a, b)
end
function LinearAlgebra._dot_nonrecursive(
    a::AbstractArray{<:Union{AbstractPolynomial,AbstractTerm}},
    b::AbstractArray{<:Number},
)
    return SA.sum_products(a, b)
end

# A bare monomial takes its coefficient from the numeric factor. Existing
# terms and polynomials instead convert numeric factors to their coefficients.
function _sum_monomial_products(a, b)
    MA._check_same_length(a, b)
    return sum(map(term, vec(a), vec(b)))
end
function LinearAlgebra._dot_nonrecursive(
    a::AbstractArray{<:Number},
    b::AbstractArray{<:AbstractMonomialLike},
)
    return _sum_monomial_products(a, b)
end
function LinearAlgebra._dot_nonrecursive(
    a::AbstractArray{<:AbstractMonomialLike},
    b::AbstractArray{<:Number},
)
    return _sum_monomial_products(b, a)
end

function _matrix_product_array(a, ::Type)
    return a
end
function _matrix_product_array(
    a::AbstractArray{<:AbstractMonomialLike},
    ::Type{T},
) where {T}
    return map(m -> term(one(T), m), a)
end

for (A, B) in (
    (AbstractMonomialLike, AbstractMonomialLike),
    (AbstractMonomialLike, Union{Number,AbstractPolynomial,AbstractTerm}),
    (Union{Number,AbstractPolynomial,AbstractTerm}, AbstractMonomialLike),
)
    @eval function MA.operate_to!(
        output::VecOrMat{P},
        ::typeof(*),
        A::AbstractMatrix{<:$A},
        B::AbstractVecOrMat{<:$B},
        α::Number = true,
    ) where {P<:AbstractPolynomial}
        # Bare monomials acquire the numeric product's coefficient type;
        # algebra factors keep their own coefficient types.
        a = _matrix_product_array(
            A,
            eltype(B) <: Number ? coefficient_type(P) : Int,
        )
        b = _matrix_product_array(
            B,
            eltype(A) <: Number ? coefficient_type(P) : Int,
        )
        return MA.operate_to!(output, *, a, b, α)
    end
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
Base.:^(x::AbstractPolynomialLike, p::Integer) = Base.power_by_squaring(x, p)
# ^(::SA.Term, ::Integer) is defined in StarAlgebras

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

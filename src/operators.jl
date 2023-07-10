# We reverse the order of comparisons here so that the result
# of x < y is equal to the result of Monomial(x) < Monomial(y)
Base.@pure function Base.isless(v1::AbstractVariable, v2::AbstractVariable)
    return name(v1) > name(v2)
end
function Base.isless(m1::AbstractTermLike, m2::AbstractTermLike)
    return isless(promote(m1, m2)...)
end

# Implement this to make coefficients be compared with terms.
function isless_coefficient(a::Real, b::Real)
    return a < b
end
function isless_coefficient(a::Number, b::Number)
    return abs(a) < abs(b)
end
# By default, coefficients are not comparable so `a` is not strictly
# less than `b`, they are considered sort of equal.
isless_coefficient(a, b) = false

function Base.isless(t1::AbstractTerm, t2::AbstractTerm)
    if monomial(t1) < monomial(t2)
        return true
    elseif monomial(t1) == monomial(t2)
        return isless_coefficient(coefficient(t1), coefficient(t2))
    else
        return false
    end
end

# promoting multiplication is not a good idea
# For example a polynomial of Float64 * a polynomial of JuMP affine expression
# is a polynomial of JuMP affine expression but if we promote it would be a
# polynomial of quadratic expression
for op in [:+, :-, :(==)]
    @eval Base.$op(p1::_APL, p2::_APL) = $op(promote(p1, p2)...)
end
# Promotion between `I` and `1` is `Any`.
# Promotion between `I` and `2I` is `UniformScaling`.
for op in [:+, :-]
    @eval function Base.$op(p1::_APL, p2::_APL{<:LinearAlgebra.UniformScaling})
        return $op(p1, map_coefficients(J -> J.λ, p2, nonzero = true))
    end
    @eval function Base.$op(p1::_APL{<:LinearAlgebra.UniformScaling}, p2::_APL)
        return $op(map_coefficients(J -> J.λ, p1, nonzero = true), p2)
    end
    @eval function Base.$op(
        p1::_APL{<:LinearAlgebra.UniformScaling},
        p2::_APL{<:LinearAlgebra.UniformScaling},
    )
        return $op(map_coefficients(J -> J.λ, p1, nonzero = true), p2)
    end
end
function Base.isapprox(p1::_APL, p2::_APL; kwargs...)
    return isapprox(promote(p1, p2)...; kwargs...)
end

# @eval $op(p::_APL, α) = $op(promote(p, α)...) would be less efficient
for (op, fun) in [
    (:+, :right_constant_plus),
    (:-, :right_constant_minus),
    (:*, :right_constant_mult),
    (:(==), :right_constant_eq),
]
    @eval Base.$op(p::_APL, α) = $fun(p, α)
end
for (op, fun) in [
    (:+, :left_constant_plus),
    (:-, :left_constant_minus),
    (:*, :left_constant_mult),
    (:(==), :left_constant_eq),
]
    @eval Base.$op(α, p::_APL) = $fun(α, p)
end
## Fix ambiguity between above methods and methods in MA
Base.:+(::MA.Zero, p::_APL) = MA.copy_if_mutable(p)
Base.:+(p::_APL, ::MA.Zero) = MA.copy_if_mutable(p)
Base.:-(::MA.Zero, p::_APL) = MA.operate(-, p)
Base.:-(p::_APL, ::MA.Zero) = MA.copy_if_mutable(p)

# Special case AbstractArrays of _APLs
# We add these instead of relying on the broadcasting API since the above method definitions are very wide.
# In particular, there is support for Matrices as coefficents. In order to avoid issues like #104 we therefore
# explicitly define this (instead of implictly getting unexpected results).
for op in [:+, :-]
    @eval Base.$op(p::_APL, A::AbstractArray{<:_APL}) = map(f -> $op(p, f), A)
    @eval Base.$op(A::AbstractArray{<:_APL}, p::_APL) = map(f -> $op(f, p), A)
end
Base.:*(p::_APL, A::AbstractArray) = map(f -> p * f, A)
Base.:*(A::AbstractArray, p::_APL) = map(f -> f * p, A)
Base.:/(A::AbstractArray, p::_APL) = map(f -> f / p, A)

right_constant_function(::typeof(+)) = right_constant_plus
right_constant_function(::typeof(-)) = right_constant_minus
right_constant_function(::typeof(*)) = right_constant_mult
function MA.operate!(op::Union{typeof(+),typeof(-),typeof(*)}, p::_APL, α)
    return MA.operate!(right_constant_function(op), p, α)
end

MA.operate!(op::typeof(*), α, p::_APL) = MA.operate!(left_constant_mult, α, p)
MA.operate!(op::typeof(*), p::_APL, α) = MA.operate!(right_constant_mult, p, α)
MA.operate!(op::typeof(/), p::_APL, α) = map_coefficients!(Base.Fix2(op, α), p)
function MA.operate_to!(output::AbstractPolynomial, op::typeof(*), α, p::_APL)
    return MA.operate_to!(output, left_constant_mult, α, p)
end
function MA.operate_to!(output::AbstractPolynomial, op::typeof(*), p::_APL, α)
    return MA.operate_to!(output, right_constant_mult, p, α)
end
function MA.operate_to!(output::_APL, op::typeof(/), p::_APL, α)
    return map_coefficients_to!(output, Base.Fix2(op, α), p)
end

function polynomial_merge!(
    n1::Int,
    n2::Int,
    get1::F1,
    get2::F2,
    set::F3,
    push::F4,
    compare_monomials::F5,
    combine::F6,
    keep::F7,
    resize::F8,
) where {F1,F2,F3,F4,F5,F6,F7,F8}
    buffer = nothing
    i = j = k = 1
    # Invariant:
    # The terms p[0] -> p[k-1] are sorted and are smaller than the remaining terms.
    # The terms p[k] -> p[i-1] are garbage.
    # The terms p[i] -> p[end] are sorted and still need to be added.
    # The terms q[j] -> p[end] are sorted and still need to be added.
    # If `buffer` is not empty:
    #   The terms in `buffer` are sorted and still need to be added.
    #   Moreover, they are smaller than the terms p[i] -> p[end].
    while i <= n1 && j <= n2
        @assert buffer === nothing || isempty(buffer)
        comp = compare_monomials(i, j)
        if comp < 0
            if k == i
                t0 = get1(i)
                if buffer === nothing
                    buffer = DataStructures.Queue{typeof(t0)}()
                end
                DataStructures.enqueue!(buffer, t0)
                i += 1
            end
            set(k, get2(j))
            j += 1
            k += 1
        elseif iszero(comp)
            combine(i, j)
            if keep(i)
                if k != i
                    @assert k < i
                    set(k, get1(i))
                end
                k += 1
            end
            i += 1
            j += 1
        else
            if k != i
                set(k, get1(i))
            end
            i += 1
            k += 1
        end
        while buffer !== nothing && !isempty(buffer) && j <= n2
            @assert i == k
            t = first(buffer)
            comp = compare_monomials(t, j)
            if comp <= 0
                if comp < 0
                    t = get2(j)
                else
                    t = combine(t, j)
                end
                j += 1
            end
            if comp >= 0
                DataStructures.dequeue!(buffer)
            end
            # if `comp` is zero, we called `combine` so `t`
            # might not be kept. If `comp` is not zero, we
            # skip the `keep` call that might be costly.
            if iszero(comp) && !keep(t)
                continue
            end
            if k <= n1
                DataStructures.enqueue!(buffer, get1(i))
                set(k, t)
            else
                push(t)
                n1 += 1
            end
            i += 1
            k += 1
        end
    end
    if buffer !== nothing && !isempty(buffer)
        @assert j == n2 + 1
        @assert i == k
        n = length(buffer)
        resize(n1 + n)
        for k in n1:-1:i
            set(k + n, get1(k))
        end
        for k in i:(i+n-1)
            set(k, DataStructures.dequeue!(buffer))
        end
        n1 += n
    else
        len = (k - 1) + (n2 - (j - 1)) + (n1 - (i - 1))
        if n1 < len
            resize(len)
        end
        while j <= n2
            set(k, get2(j))
            j += 1
            k += 1
        end
        @assert j == n2 + 1
        while i <= n1
            set(k, get1(i))
            i += 1
            k += 1
        end
        if len < n1
            resize(len)
        end
        @assert i == n1 + 1
        @assert k == len + 1
    end
    return
end

#MA.operate!(op::Union{typeof(+), typeof(-)}, p::AbstractPolynomial, q::AbstractPolynomial) = MA.operate_to!(p, op, p, q)
function MA.operate!(
    op::Union{typeof(+),typeof(-)},
    p::AbstractPolynomial,
    q::AbstractPolynomialLike,
)
    return MA.operate!(op, p, polynomial(q))
end

function mul_to_terms!(ts::Vector{<:AbstractTerm}, p1::_APL, p2::_APL)
    for t1 in terms(p1)
        for t2 in terms(p2)
            push!(ts, t1 * t2)
        end
    end
    return ts
end
function Base.:*(p::AbstractPolynomial, q::AbstractPolynomial)
    return polynomial!(
        mul_to_terms!(
            MA.promote_operation(*, term_type(p), term_type(q))[],
            p,
            q,
        ),
    )
end

Base.isapprox(p::_APL, α; kwargs...) = isapprox(promote(p, α)...; kwargs...)
Base.isapprox(α, p::_APL; kwargs...) = isapprox(promote(p, α)...; kwargs...)

# `MA.operate(-, p)` redirects to `-p` as it assumes that `-p` can be modified
# through the MA API without modifying `p`. We should either copy the monomial
# here or implement a `MA.operate(-, p)` that copies it. We choose the first
# option.
Base.:-(m::AbstractMonomialLike) = _term(-1, MA.copy_if_mutable(m))
Base.:-(t::AbstractTermLike) = _term(MA.operate(-, coefficient(t)), monomial(t))
Base.:-(p::_APL) = map_coefficients(-, p)
Base.:+(p::Union{_APL,RationalPoly}) = p
Base.:*(p::Union{_APL,RationalPoly}) = p

# Avoid adding a zero constant that might artificially increase the Newton polytope
# Need to add polynomial conversion for type stability
function right_constant_plus(p::_APL{S}, α::T) where {S,T}
    return iszero(α) ? polynomial(p, MA.promote_operation(+, S, T)) :
           p + constant_term(α, p)
end
function left_constant_plus(α::S, p::_APL{T}) where {S,T}
    return iszero(α) ? polynomial(p, MA.promote_operation(+, S, T)) :
           constant_term(α, p) + p
end
function MA.operate!(::typeof(right_constant_plus), p::_APL, α)
    if !iszero(α)
        MA.operate!(+, p, constant_term(α, p))
    end
    return p
end
function right_constant_minus(p::_APL{S}, α::T) where {S,T}
    return iszero(α) ? polynomial(p, MA.promote_operation(-, S, T)) :
           p - constant_term(α, p)
end
function left_constant_minus(α::S, p::_APL{T}) where {S,T}
    return iszero(α) ? polynomial(-p, MA.promote_operation(-, S, T)) :
           constant_term(α, p) - p
end
function MA.operate!(::typeof(right_constant_minus), p::_APL, α)
    if !iszero(α)
        MA.operate!(-, p, constant_term(α, p))
    end
    return p
end

# Coefficients and variables commute
left_constant_mult(α, v::AbstractMonomial) = term_type(v, typeof(α))(α, v)
left_constant_mult(α, v::AbstractVariable) = left_constant_mult(α, monomial(v)) # TODO linear term
right_constant_mult(m::AbstractMonomialLike, α) = left_constant_mult(α, m)

function left_constant_mult(α, p::AbstractPolynomialLike)
    return map_coefficients(Base.Fix1(*, α), p)
end
function right_constant_mult(p::AbstractPolynomialLike, α)
    return map_coefficients(Base.Fix2(*, α), p)
end

function MA.operate_to!(output, ::typeof(left_constant_mult), α, p::_APL)
    return map_coefficients_to!(output, Base.Fix1(*, α), p)
end
function MA.operate_to!(output, ::typeof(right_constant_mult), p::_APL, α)
    return map_coefficients_to!(output, Base.Fix2(*, α), p)
end
function MA.operate!(::typeof(left_constant_mult), α, p::_APL)
    return map_coefficients!(Base.Fix1(*, α), p)
end
function MA.operate!(::typeof(right_constant_mult), p::_APL, α)
    return map_coefficients!(Base.Fix2(MA.mul!!, α), p)
end

function MA.operate_to!(
    output::AbstractMonomial,
    ::typeof(*),
    m1::AbstractMonomialLike,
    m2::AbstractMonomialLike,
)
    return map_exponents_to!(output, +, m1, m2)
end
function MA.operate!(
    ::typeof(*),
    m1::AbstractMonomial,
    m2::AbstractMonomialLike,
)
    return map_exponents!(+, m1, m2)
end
function Base.:*(m1::AbstractMonomialLike, m2::AbstractMonomialLike)
    return map_exponents(+, m1, m2)
end
#Base.:*(m1::AbstractMonomialLike, m2::AbstractMonomialLike) = *(monomial(m1), monomial(m2))

function Base.:*(m::AbstractMonomialLike, t::AbstractTermLike)
    return left_constant_mult(coefficient(t), m * monomial(t))
end
function Base.:*(t::AbstractTermLike, m::AbstractMonomialLike)
    return left_constant_mult(coefficient(t), monomial(t) * m)
end
function Base.:*(t1::AbstractTermLike, t2::AbstractTermLike)
    return left_constant_mult(
        coefficient(t1) * coefficient(t2),
        monomial(t1) * monomial(t2),
    )
end

function MA.operate!(::typeof(*), p::_APL, t::AbstractMonomialLike)
    return map_exponents!(+, p, t)
end
Base.:*(p::_APL, t::AbstractMonomialLike) = map_exponents(+, p, t)
Base.:*(t::AbstractTermLike, p::_APL) = polynomial!(map(te -> t * te, terms(p)))
Base.:*(p::_APL, t::AbstractTermLike) = polynomial!(map(te -> te * t, terms(p)))
Base.:*(p::_APL, q::_APL) = polynomial(p) * polynomial(q)

# guaranteed that monomial(t1) > monomial(t2)
function _polynomial_2terms(
    t1::TT,
    t2::TT,
    ::Type{T},
) where {TT<:AbstractTerm,T}
    if iszero(t1)
        polynomial(t2, T)
    elseif iszero(t2)
        polynomial(t1, T)
    else
        # not `polynomial!` because we `t1` and `t2` cannot be modified
        polynomial(term_type(TT, T)[t1, t2], SortedUniqState())
    end
end

_term(α, mono) = term(α, MA.copy_if_mutable(mono))

for op in [:+, :-]
    @eval begin
        function Base.$op(t1::AbstractTermLike, t2::AbstractTermLike)
            return $op(term(t1), term(t2))
        end
        function Base.$op(t1::AbstractTerm, t2::AbstractTerm)
            return $op(_promote_terms(t1, t2)...)
        end
        function Base.$op(t1::TT, t2::TT) where {T,TT<:AbstractTerm{T}}
            S = MA.promote_operation($op, T, T)
            # t1 > t2 would compare the coefficient in case the monomials are equal
            # and it will throw a MethodError in case the coefficients are not comparable
            if monomial(t1) == monomial(t2)
                return polynomial(
                    _term($op(coefficient(t1), coefficient(t2)), monomial(t1)),
                    S,
                )
            elseif monomial(t1) < monomial(t2)
                return _polynomial_2terms(t1, $op(t2), S)
            else
                return _polynomial_2terms($op(t2), t1, S)
            end
        end
    end
end
_promote_terms(t1, t2) = promote(t1, t2)
# Promotion between `I` and `1` is `Any`.
function _promote_terms(
    t1::AbstractTerm,
    t2::AbstractTerm{<:LinearAlgebra.UniformScaling},
)
    return _promote_terms(t1, coefficient(t2).λ * monomial(t2))
end
function _promote_terms(
    t1::AbstractTerm{<:LinearAlgebra.UniformScaling},
    t2::AbstractTerm,
)
    return _promote_terms(coefficient(t1).λ * monomial(t1), t2)
end
# Promotion between `I` and `2I` is `UniformScaling`, not `UniformScaling{Int}`.
function _promote_terms(
    t1::AbstractTerm{LinearAlgebra.UniformScaling{S}},
    t2::AbstractTerm{LinearAlgebra.UniformScaling{T}},
) where {S<:Number,T<:Number}
    U = LinearAlgebra.UniformScaling{promote_type(S, T)}
    return _promote_terms(
        MA.scaling_convert(U, coefficient(t1)) * monomial(t1),
        MA.scaling_convert(U, coefficient(t2)) * monomial(t2),
    )
end
function _promote_terms(
    t1::AbstractTerm{LinearAlgebra.UniformScaling{T}},
    t2::AbstractTerm{LinearAlgebra.UniformScaling{T}},
) where {T<:Number}
    return promote(t1, t2)
end

LinearAlgebra.adjoint(v::AbstractVariable) = conj(v)
LinearAlgebra.adjoint(m::AbstractMonomial) = conj(m)
function LinearAlgebra.adjoint(t::AbstractTerm)
    return _term(adjoint(coefficient(t)), adjoint(monomial(t)))
end
function LinearAlgebra.adjoint(p::AbstractPolynomialLike)
    return polynomial!(adjoint.(terms(p)))
end
function LinearAlgebra.adjoint(r::RationalPoly)
    return adjoint(numerator(r)) / adjoint(denominator(r))
end
LinearAlgebra.hermitian_type(::Type{T}) where {T<:AbstractPolynomialLike} = T
function LinearAlgebra.hermitian(v::AbstractPolynomialLike, ::Symbol)
    iscomplex(v) && error(
        "Complex-valued polynomials cannot be interpreted as hermitian scalars",
    )
    return v
end

LinearAlgebra.transpose(v::AbstractVariable) = v
LinearAlgebra.transpose(m::AbstractMonomial) = m
function LinearAlgebra.transpose(t::AbstractTerm)
    return _term(LinearAlgebra.transpose(coefficient(t)), monomial(t))
end
function LinearAlgebra.transpose(p::AbstractPolynomialLike)
    return polynomial(map(LinearAlgebra.transpose, terms(p)))
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
LinearAlgebra.dot(x, p::AbstractPolynomialLike) = x' * p
LinearAlgebra.dot(p::AbstractPolynomialLike, x) = p' * x

LinearAlgebra.symmetric_type(PT::Type{<:_APL}) = PT
LinearAlgebra.symmetric(p::_APL, ::Symbol) = p

# Amazingly, this works! Thanks, StaticArrays.jl!
"""
Convert a tuple of variables into a static vector to allow array-like usage.
The element type of the vector will be Monomial{vars, length(vars)}.
"""
Base.vec(vars::Tuple{Vararg{AbstractVariable}}) = [vars...]
# vec(vars::Tuple{Vararg{<:TypedVariable}}) = SVector(vars)

# https://github.com/JuliaLang/julia/pull/23332
Base.:^(x::AbstractPolynomialLike, p::Integer) = Base.power_by_squaring(x, p)

function MA.operate_to!(
    output::AbstractPolynomial,
    op::MA.AddSubMul,
    x,
    args::Vararg{Any,N},
) where {N}
    return MA.operate_to!(output, MA.add_sub_op(op), x, *(args...))
end
function MA.operate!(
    op::MA.AddSubMul,
    x::AbstractPolynomial,
    y,
    z,
    args::Vararg{Any,N},
) where {N}
    return MA.operate!(MA.add_sub_op(op), x, *(y, z, args...))
end
function MA.buffer_for(
    ::MA.AddSubMul,
    ::Type{<:AbstractPolynomial},
    args::Vararg{Type,N},
) where {N}
    return zero(MA.promote_operation(*, args...))
end
function MA.buffered_operate_to!(
    buffer::AbstractPolynomial,
    output::AbstractPolynomial,
    op::MA.AddSubMul,
    x::AbstractPolynomial,
    y,
    z,
    args::Vararg{Any,N},
) where {N}
    product = MA.operate_to!!(buffer, *, y, z, args...)
    return MA.operate_to!(output, MA.add_sub_op(op), x, product)
end
function MA.buffered_operate!(
    buffer,
    op::MA.AddSubMul,
    x::AbstractPolynomial,
    y,
    z,
    args::Vararg{Any,N},
) where {N}
    product = MA.operate_to!!(buffer, *, y, z, args...)
    return MA.operate!(MA.add_sub_op(op), x, product)
end

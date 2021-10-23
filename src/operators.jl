# We reverse the order of comparisons here so that the result
# of x < y is equal to the result of Monomial(x) < Monomial(y)
Base.@pure Base.isless(v1::AbstractVariable, v2::AbstractVariable) = name(v1) > name(v2)
Base.isless(m1::AbstractTermLike, m2::AbstractTermLike) = isless(promote(m1, m2)...)

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
    @eval Base.$op(p1::APL, p2::APL) = $op(promote(p1, p2)...)
end
# Promotion between `I` and `1` is `Any`.
# Promotion between `I` and `2I` is `UniformScaling`.
for op in [:+, :-]
    @eval Base.$op(p1::APL, p2::APL{<:LinearAlgebra.UniformScaling}) = $op(p1, mapcoefficientsnz(J -> J.λ, p2))
    @eval Base.$op(p1::APL{<:LinearAlgebra.UniformScaling}, p2::APL) = $op(mapcoefficientsnz(J -> J.λ, p1), p2)
    @eval Base.$op(p1::APL{<:LinearAlgebra.UniformScaling}, p2::APL{<:LinearAlgebra.UniformScaling}) = $op(mapcoefficientsnz(J -> J.λ, p1), p2)
end
Base.isapprox(p1::APL, p2::APL; kwargs...) = isapprox(promote(p1, p2)...; kwargs...)

# @eval $op(p::APL, α) = $op(promote(p, α)...) would be less efficient
for (op, fun) in [(:+, :plusconstant), (:-, :minusconstant), (:*, :multconstant), (:(==), :eqconstant)]
    @eval Base.$op(p::APL, α) = $fun(p, α)
    @eval Base.$op(α, p::APL) = $fun(α, p)
end
# Fix ambiguity between above methods and methods in MA
Base.:+(::MA.Zero, p::APL) = MA.copy_if_mutable(p)
Base.:+(p::APL, ::MA.Zero) = MA.copy_if_mutable(p)
Base.:-(::MA.Zero, p::APL) = MA.operate(-, p)
Base.:-(p::APL, ::MA.Zero) = MA.copy_if_mutable(p)

# Special case AbstractArrays of APLs
# We add these instead of relying on the broadcasting API since the above method definitions are very wide.
# In particular, there is support for Matrices as coefficents. In order to avoid issues like #104 we therefore
# explicitly define this (instead of implictly getting unexpected results).
for op in [:+, :-]
    @eval Base.$op(p::APL, A::AbstractArray{<:APL}) = map(f -> $op(p, f), A)
    @eval Base.$op(A::AbstractArray{<:APL}, p::APL) = map(f -> $op(f, p), A)
end
Base.:*(p::APL, A::AbstractArray) = map(f -> p * f, A)
Base.:*(A::AbstractArray, p::APL) = map(f -> f * p, A)
Base.:/(A::AbstractArray, p::APL) = map(f -> f / p, A)

constant_function(::typeof(+)) = plusconstant
constant_function(::typeof(-)) = minusconstant
MA.operate!(op::Union{typeof(+), typeof(-)}, p::APL, α) = MA.operate!(constant_function(op), p, α)
MA.operate_to!(output::AbstractPolynomial, op::typeof(*), α, p::APL) = MA.operate_to!(output, multconstant, α, p)
MA.operate_to!(output::AbstractPolynomial, op::typeof(*), p::APL, α) = MA.operate_to!(output, multconstant, p, α)

function polynomial_merge!(
    n1::Int, n2::Int, get1::Function, get2::Function,
    set::Function, push::Function, compare_monomials::Function,
    combine::Function, keep::Function, resize::Function)
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
        if comp > 0
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
            if comp >= 0
                if comp > 0
                    t = get2(j)
                else
                    t = combine(t, j)
                end
                j += 1
            end
            if comp <= 0
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
        for k in i:(i + n - 1)
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
MA.operate!(op::Union{typeof(+), typeof(-)}, p::AbstractPolynomial, q::AbstractPolynomialLike) = MA.operate!(op, p, polynomial(q))

function mul_to_terms!(ts::Vector{<:AbstractTerm}, p1::APL, p2::APL)
    for t1 in terms(p1)
        for t2 in terms(p2)
            push!(ts, t1 * t2)
        end
    end
    return ts
end
function Base.:*(p::AbstractPolynomial, q::AbstractPolynomial)
    polynomial!(mul_to_terms!(MA.promote_operation(*, termtype(p), termtype(q))[], p, q))
end

Base.isapprox(p::APL, α; kwargs...) = isapprox(promote(p, α)...; kwargs...)
Base.isapprox(α, p::APL; kwargs...) = isapprox(promote(p, α)...; kwargs...)

# `MA.operate(-, p)` redirects to `-p` as it assumes that `-p` can be modified
# through the MA API without modifying `p`. We should either copy the monomial
# here or implement a `MA.operate(-, p)` that copies it. We choose the first
# option.
Base.:-(m::AbstractMonomialLike) = _term(-1, MA.copy_if_mutable(m))
Base.:-(t::AbstractTermLike) = _term(MA.operate(-, coefficient(t)), monomial(t))
Base.:-(p::APL) = polynomial!((-).(terms(p)))
Base.:+(p::Union{APL, RationalPoly}) = p
Base.:*(p::Union{APL, RationalPoly}) = p

# Avoid adding a zero constant that might artificially increase the Newton polytope
# Need to add polynomial conversion for type stability
plusconstant(p::APL{S}, α::T)  where {S, T} = iszero(α) ? polynomial( p, MA.promote_operation(+, S, T)) : p + constantterm(α, p)
plusconstant(α::S, p::APL{T})  where {S, T} = iszero(α) ? polynomial( p, MA.promote_operation(+, S, T)) : constantterm(α, p) + p
function MA.operate!(::typeof(plusconstant), p::APL, α)
    if !iszero(α)
        MA.operate!(+, p, constantterm(α, p))
    end
    return p
end
minusconstant(p::APL{S}, α::T) where {S, T} = iszero(α) ? polynomial( p, MA.promote_operation(-, S, T)) : p - constantterm(α, p)
minusconstant(α::S, p::APL{T}) where {S, T} = iszero(α) ? polynomial(-p, MA.promote_operation(-, S, T)) : constantterm(α, p) - p
function MA.operate!(::typeof(minusconstant), p::APL, α)
    if !iszero(α)
        MA.operate!(-, p, constantterm(α, p))
    end
    return p
end

# Coefficients and variables commute
multconstant(α, v::AbstractVariable) = multconstant(α, monomial(v)) # TODO linear term
multconstant(m::AbstractMonomialLike, α) = multconstant(α, m)

_multconstant(α, f, t::AbstractTermLike) = mapcoefficientsnz(f, t)
function _multconstant(α::T, f, p::AbstractPolynomial{S}) where {S, T}
    if iszero(α)
        zero(polynomialtype(p, MA.promote_operation(*, T, S)))
    else
        mapcoefficientsnz(f, p)
    end
end
_multconstant(α, f, p::AbstractPolynomialLike) = _multconstant(α, f, polynomial(p))

multconstant(α, p::AbstractPolynomialLike) = _multconstant(α, β -> α*β, p)
multconstant(p::AbstractPolynomialLike, α) = _multconstant(α, β -> β*α, p)

function mapcoefficientsnz_to! end

function _multconstant_to!(output, α, f, p)
    if iszero(α)
        MA.operate!(zero, output)
    else
        mapcoefficientsnz_to!(output, f, p)
    end
end
function MA.operate_to!(output, ::typeof(multconstant), p::APL, α)
    _multconstant_to!(output, α, β -> β*α, p)
end
function MA.operate_to!(output, ::typeof(multconstant), α, p::APL)
    _multconstant_to!(output, α, β -> α*β, p)
end

MA.operate_to!(output::AbstractMonomial, ::typeof(*), m1::AbstractMonomialLike, m2::AbstractMonomialLike) = mapexponents_to!(output, +, m1, m2)
MA.operate!(::typeof(*), m1::AbstractMonomial, m2::AbstractMonomialLike) = mapexponents!(+, m1, m2)
Base.:*(m1::AbstractMonomialLike, m2::AbstractMonomialLike) = mapexponents(+, m1, m2)
#Base.:*(m1::AbstractMonomialLike, m2::AbstractMonomialLike) = *(monomial(m1), monomial(m2))

Base.:*(m::AbstractMonomialLike, t::AbstractTermLike) = term(coefficient(t), m * monomial(t))
Base.:*(t::AbstractTermLike, m::AbstractMonomialLike) = term(coefficient(t), monomial(t) * m)
Base.:*(t1::AbstractTermLike, t2::AbstractTermLike) = term(coefficient(t1) * coefficient(t2), monomial(t1) * monomial(t2))

Base.:*(t::AbstractTermLike, p::APL) = polynomial!(map(te -> t * te, terms(p)))
Base.:*(p::APL, t::AbstractTermLike) = polynomial!(map(te -> te * t, terms(p)))
Base.:*(p::APL, q::APL) = polynomial(p) * polynomial(q)

# guaranteed that monomial(t1) > monomial(t2)
function _polynomial_2terms(t1::TT, t2::TT, ::Type{T}) where {TT<:AbstractTerm, T}
    if iszero(t1)
        polynomial(t2, T)
    elseif iszero(t2)
        polynomial(t1, T)
    else
        # not `polynomial!` because we `t1` and `t2` cannot be modified
        polynomial(termtype(TT, T)[t1, t2], SortedUniqState())
    end
end

_term(α, mono) = term(α, MA.copy_if_mutable(mono))

for op in [:+, :-]
    @eval begin
        Base.$op(t1::AbstractTermLike, t2::AbstractTermLike) = $op(term(t1), term(t2))
        Base.$op(t1::AbstractTerm, t2::AbstractTerm) = $op(_promote_terms(t1, t2)...)
        function Base.$op(t1::TT, t2::TT) where {T, TT <: AbstractTerm{T}}
            S = MA.promote_operation($op, T, T)
            # t1 > t2 would compare the coefficient in case the monomials are equal
            # and it will throw a MethodError in case the coefficients are not comparable
            if monomial(t1) == monomial(t2)
                polynomial(_term($op(coefficient(t1), coefficient(t2)), monomial(t1)), S)
            elseif monomial(t1) > monomial(t2)
                ts = _polynomial_2terms(t1, $op(t2), S)
            else
                ts = _polynomial_2terms($op(t2), t1, S)
            end
        end
    end
end
_promote_terms(t1, t2) = promote(t1, t2)
# Promotion between `I` and `1` is `Any`.
_promote_terms(t1::AbstractTerm, t2::AbstractTerm{<:LinearAlgebra.UniformScaling}) = _promote_terms(t1, coefficient(t2).λ * monomial(t2))
_promote_terms(t1::AbstractTerm{<:LinearAlgebra.UniformScaling}, t2::AbstractTerm) = _promote_terms(coefficient(t1).λ * monomial(t1), t2)
# Promotion between `I` and `2I` is `UniformScaling`, not `UniformScaling{Int}`.
function _promote_terms(t1::AbstractTerm{LinearAlgebra.UniformScaling{S}}, t2::AbstractTerm{LinearAlgebra.UniformScaling{T}}) where {S<:Number, T<:Number}
    U = LinearAlgebra.UniformScaling{promote_type(S, T)}
    _promote_terms(MA.scaling_convert(U, coefficient(t1)) * monomial(t1), MA.scaling_convert(U, coefficient(t2)) * monomial(t2))
end
function _promote_terms(t1::AbstractTerm{LinearAlgebra.UniformScaling{T}}, t2::AbstractTerm{LinearAlgebra.UniformScaling{T}}) where T<:Number
    return promote(t1, t2)
end

LinearAlgebra.adjoint(v::AbstractVariable) = v
LinearAlgebra.adjoint(m::AbstractMonomial) = m
LinearAlgebra.adjoint(t::AbstractTerm) = _term(LinearAlgebra.adjoint(coefficient(t)), monomial(t))
LinearAlgebra.adjoint(p::AbstractPolynomialLike) = polynomial(map(LinearAlgebra.adjoint, terms(p)))
LinearAlgebra.adjoint(r::RationalPoly) = adjoint(numerator(r)) / adjoint(denominator(r))

LinearAlgebra.transpose(v::AbstractVariable) = v
LinearAlgebra.transpose(m::AbstractMonomial) = m
LinearAlgebra.transpose(t::AbstractTerm) = _term(LinearAlgebra.transpose(coefficient(t)), monomial(t))
LinearAlgebra.transpose(p::AbstractPolynomialLike) = polynomial(map(LinearAlgebra.transpose, terms(p)))
LinearAlgebra.transpose(r::RationalPoly) = transpose(numerator(r)) / transpose(denominator(r))

LinearAlgebra.dot(p1::AbstractPolynomialLike, p2::AbstractPolynomialLike) = p1' * p2
LinearAlgebra.dot(x, p::AbstractPolynomialLike) = x' * p
LinearAlgebra.dot(p::AbstractPolynomialLike, x) = p' * x

LinearAlgebra.symmetric_type(PT::Type{<:APL}) = PT
LinearAlgebra.symmetric(p::APL, ::Symbol) = p

# Amazingly, this works! Thanks, StaticArrays.jl!
"""
Convert a tuple of variables into a static vector to allow array-like usage.
The element type of the vector will be Monomial{vars, length(vars)}.
"""
Base.vec(vars::Tuple{Vararg{AbstractVariable}}) = [vars...]
# vec(vars::Tuple{Vararg{<:TypedVariable}}) = SVector(vars)

# https://github.com/JuliaLang/julia/pull/23332
Base.:^(x::AbstractPolynomialLike, p::Integer) = Base.power_by_squaring(x, p)

function MA.operate_to!(output::AbstractPolynomial, op::MA.AddSubMul, x, args::Vararg{Any, N}) where N
    return MA.operate_to!(output, MA.add_sub_op(op), x, *(args...))
end
function MA.operate!(op::MA.AddSubMul, x::AbstractPolynomial, y, z, args::Vararg{Any, N}) where N
    return MA.operate!(MA.add_sub_op(op), x, *(y, z, args...))
end
MA.buffer_for(::MA.AddSubMul, ::Type{<:AbstractPolynomial}, args::Vararg{Type, N}) where {N} = zero(MA.promote_operation(*, args...))
function MA.buffered_operate_to!(buffer::AbstractPolynomial, output::AbstractPolynomial, op::MA.AddSubMul, x::AbstractPolynomial, y, z, args::Vararg{Any, N}) where N
    product = MA.operate_to!!(buffer, *, y, z, args...)
    return MA.operate_to!(output, MA.add_sub_op(op), x, product)
end
function MA.buffered_operate!(buffer::AbstractPolynomial, op::MA.AddSubMul, x::AbstractPolynomial, y, z, args::Vararg{Any, N}) where N
    product = MA.operate_to!!(buffer, *, y, z, args...)
    return MA.operate!(MA.add_sub_op(op), x, product)
end

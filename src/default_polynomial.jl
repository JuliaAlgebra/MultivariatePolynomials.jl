struct Polynomial{CoeffType,T<:AbstractTerm{CoeffType},V<:AbstractVector{T}} <: AbstractPolynomial{CoeffType}
    terms::V

    Polynomial{C, T, V}(terms::AbstractVector{T}) where {C, T, V} = new{C, T, V}(terms)
end

function coefficients(p::Polynomial)
    return LazyMap{coefficienttype(p)}(coefficient, terms(p))
end
function monomials(p::Polynomial)
    return LazyMap{monomial_type(p)}(monomial, terms(p))
end

const VectorPolynomial{C,T} = Polynomial{C,T,Vector{T}}

term_type(::Type{<:Polynomial{C,T}}) where {C,T} = T
terms(p::Polynomial) = p.terms
constant_monomial(::Union{Polynomial{C,TT},Type{Polynomial{C,TT}}}) where {C,TT} = constant_monomial(TT)
Base.convert(::Type{Polynomial{C,TT,V}}, p::Polynomial{C,TT,V}) where {C,TT,V} = p
function Base.convert(PT::Type{Polynomial{C,TT,V}}, p::AbstractPolynomialLike) where {C,TT,V}
    return PT(convert(V, terms(p)))
end
function Base.convert(::Type{Polynomial{C}}, p::AbstractPolynomialLike) where {C}
    TT = term_type(p, C)
    return convert(Polynomial{C,TT,Vector{TT}}, p)
end
function Base.convert(::Type{Polynomial}, p::AbstractPolynomialLike)
    return convert(Polynomial{coefficienttype(p)}, p)
end


_change_eltype(::Type{<:Vector}, ::Type{T}) where T = Vector{T}
function polynomial_type(::Type{Polynomial{C, T, V}}, ::Type{NewC}) where {C, T, V, NewC}
    NewT = term_type(T, NewC)
    Polynomial{NewC, NewT, _change_eltype(V, NewT)}
end

function Base.promote_rule(::Type{Polynomial}, ::Type{<:AbstractPolynomialLike})
    return Polynomial
end
function Base.promote_rule(::Type{<:AbstractPolynomialLike}, ::Type{Polynomial})
    return Polynomial
end
function promote_rule_constant(::Type{T}, ::Type{Polynomial}) where T
    return Any
    # `convert(Polynomial{T}, ::T)` cannot be implemented as we don't know the type of the term
    #return Polynomial{T}
end

(p::Polynomial)(s...) = substitute(Eval(), p, s)

function Base.one(::Type{P}) where {C,TT,P<:Polynomial{C,TT}}
    return convert(P, one(TT))
end
Base.one(p::Polynomial) = one(typeof(p))

Base.zero(::Type{Polynomial{C,T,A}}) where {C,T,A} = Polynomial{C,T,A}(A())
Base.zero(t::Polynomial) = zero(typeof(t))

join_terms(terms1::AbstractArray{<:Term}, terms2::AbstractArray{<:Term}) = Sequences.merge_sorted(terms1, terms2, compare, combine)
function join_terms!(output::AbstractArray{<:Term}, terms1::AbstractArray{<:Term}, terms2::AbstractArray{<:Term})
    resize!(output, length(terms1) + length(terms2))
    Sequences.merge_sorted!(output, terms1, terms2, compare, combine)
end

Base.:(+)(p1::Polynomial, p2::Polynomial) = polynomial!(join_terms(terms(p1), terms(p2)), SortedUniqState())
Base.:(+)(p1::Polynomial, p2::Polynomial{<:LinearAlgebra.UniformScaling}) = p1 + map_coefficients(J -> J.λ, p2, nonzero=true)
Base.:(+)(p1::Polynomial{<:LinearAlgebra.UniformScaling}, p2::Polynomial) = map_coefficients(J -> J.λ, p1, nonzero=true) + p2
Base.:(+)(p1::Polynomial{<:LinearAlgebra.UniformScaling}, p2::Polynomial{<:LinearAlgebra.UniformScaling}) = map_coefficients(J -> J.λ, p1, nonzero=true) + p2
function MA.operate_to!(result::Polynomial, ::typeof(+), p1::Polynomial, p2::Polynomial)
    if result === p1 || result === p2
        error("Cannot call `operate_to!(output, +, p, q)` with `output` equal to `p` or `q`, call `operate!` instead.")
    end
    join_terms!(result.terms, terms(p1), terms(p2))
    result
end
function MA.operate_to!(result::Polynomial, ::typeof(*), p::Polynomial, t::AbstractTermLike)
    if iszero(t)
        MA.operate!(zero, result)
    else
        resize!(result.terms, nterms(p))
        for i in eachindex(p.terms)
            # TODO could use MA.mul_to!! for indices that were presents in `result` before the `resize!`.
            result.terms[i] = p.terms[i] * t
        end
        return result
    end
end
function MA.operate_to!(result::Polynomial, ::typeof(*), t::AbstractTermLike, p::Polynomial)
    if iszero(t)
        MA.operate!(zero, result)
    else
        resize!(result.terms, nterms(p))
        for i in eachindex(p.terms)
            # TODO could use MA.mul_to!! for indices that were presents in `result` before the `resize!`.
            result.terms[i] = t * p.terms[i]
        end
        return result
    end
end
Base.:(-)(p1::Polynomial, p2::Polynomial) = polynomial!(join_terms(terms(p1), (-).(terms(p2))))
Base.:(-)(p1::Polynomial, p2::Polynomial{<:LinearAlgebra.UniformScaling}) = p1 - map_coefficients(J -> J.λ, p2, nonzero=true)
Base.:(-)(p1::Polynomial{<:LinearAlgebra.UniformScaling}, p2::Polynomial) = map_coefficients(J -> J.λ, p1, nonzero=true) - p2
Base.:(-)(p1::Polynomial{<:LinearAlgebra.UniformScaling}, p2::Polynomial{<:LinearAlgebra.UniformScaling}) = map_coefficients(J -> J.λ, p1, nonzero=true) - p2

LinearAlgebra.adjoint(x::Polynomial) = polynomial!(adjoint.(terms(x)))

function map_coefficients(f::F, p::Polynomial; nonzero = false) where {F<:Function}
    terms = map(p.terms) do term
        map_coefficients(f, term)
    end
    if !nonzero
        filter!(!iszero, terms)
    end
    return polynomial!(terms)
end
function map_coefficients!(f::F, p::Polynomial; nonzero = false) where {F<:Function}
    for i in eachindex(p.terms)
        t = p.terms[i]
        p.terms[i] = Term(f(coefficient(t)), monomial(t))
    end
    if !nonzero
        filter!(!iszero, p.terms)
    end
    return p
end

function map_coefficients_to!(output::Polynomial, f::F, p::Polynomial; nonzero = false) where {F<:Function}
    resize!(output.terms, nterms(p))
    for i in eachindex(p.terms)
        t = p.terms[i]
        output.terms[i] = Term(f(coefficient(t)), monomial(t))
    end
    if !nonzero
        filter!(!iszero, output.terms)
    end
    return output
end
function map_coefficients_to!(output::Polynomial, f::Function, p::AbstractPolynomialLike; nonzero = false)
    return map_coefficients_to!(output, f, polynomial(p); nonzero = false)
end

# The polynomials can be mutated.
MA.mutability(::Type{<:VectorPolynomial}) = MA.IsMutable()

# Terms are not mutable under the MutableArithmetics API
function MA.mutable_copy(p::VectorPolynomial{C,TT}) where {C,TT}
    return VectorPolynomial{C,TT}([TT(MA.copy_if_mutable(coefficient(term)), monomial(term)) for term in terms(p)])
end
Base.copy(p::VectorPolynomial) = MA.mutable_copy(p)

function grlex end

function __polynomial_merge!(op::MA.AddSubMul, p::Polynomial{T,TT}, get1, set, push, resize, keep, q, t::AbstractTermLike...) where {T,TT}
    compare_monomials = let t=t, get1=get1, q=q
        (tp, j) -> begin
            if tp isa Int && j isa Int
                tp = get1(tp)
            end
            grlex(*(monomials(q)[j], monomial.(t)...), monomial(tp))
        end
    end
    get2 = let t=t, q=q
        i -> begin
            tq = terms(q)[i]
            TT(MA.scaling_convert(T, MA.operate(MA.add_sub_op(op), *(coefficient(tq), coefficient.(t)...))), *(monomial(tq), monomial.(t)...))
        end
    end
    combine = let t=t, p=p, q=q
        (i, j) -> begin
            if i isa Int
                p.terms[i] = Term(MA.operate!!(op, coefficient(p.terms[i]), coefficients(q)[j], coefficient.(t)...), monomial(p.terms[i]))
            else
                typeof(i)(MA.operate!!(op, coefficient(i), coefficients(q)[j], coefficient.(t)...), monomial(i))
            end
        end
    end
    polynomial_merge!(
        nterms(p), nterms(q), get1, get2, set, push,
        compare_monomials, combine, keep, resize
    )
    return
end

function __polynomial_merge!(op::MA.AddSubMul, p::Polynomial{T,TT}, get1, set, push, resize, keep, t::AbstractTermLike, q::Polynomial, buffer=nothing) where {T,TT}
    compare_monomials = let t=t, get1=get1, q=q
        (tp, j) -> begin
            if tp isa Int && j isa Int
                tp = get1(tp)
            end
            grlex(monomial(t) * monomials(q)[j], monomial(tp))
        end
    end
    get2 = let t=t, q=q
        i -> begin
            tq = terms(q)[i]
            TT(MA.scaling_convert(T, MA.operate(MA.add_sub_op(op), coefficient(t) * coefficient(tq))), monomial(t) * monomial(tq))
        end
    end
    combine = let t=t, p=p, q=q
        (i, j) -> begin
            if i isa Int
                p.terms[i] = Term(MA.buffered_operate!!(buffer, op, coefficient(p.terms[i]), coefficient(t), coefficients(q)[j]), monomial(p.terms[i]))
            else
                typeof(i)(MA.buffered_operate!!(buffer, op, coefficient(i), coefficient(t), coefficients(q)[j]), monomial(i))
            end
        end
    end
    polynomial_merge!(
        nterms(p), nterms(q), get1, get2, set, push,
        compare_monomials, combine, keep, resize
    )
    return
end

function __polynomial_merge!(op::Union{typeof(+), typeof(-)}, p::Polynomial{T,TT}, get1, set, push, resize, keep, q::Union{Polynomial,AbstractTermLike}) where {T,TT}
    compare_monomials = let get1=get1, q=q
        (t, j) -> begin
            if t isa Int && j isa Int
                t = get1(t)
            end
            grlex(monomials(q)[j], monomial(t))
        end
    end
    get2 = let q=q
        i -> begin
            t = terms(q)[i]
            # `operate` makes sure we make a copy of the term as it may be stored directly in `p`
            TT(MA.scaling_convert(T, MA.operate(op, coefficient(t))), monomial(t))
        end
    end
    combine = let p=p, q=q
        (i, j) -> begin
            if i isa Int
                p.terms[i] = Term(MA.operate!!(op, coefficient(p.terms[i]), coefficients(q)[j]), monomial(p.terms[i]))
            else
                typeof(i)(MA.operate!!(op, coefficient(i), coefficients(q)[j]), monomial(i))
            end
        end
    end
    polynomial_merge!(
        nterms(p), nterms(q), get1, get2, set, push,
        compare_monomials, combine, keep, resize
    )
    return
end

function _polynomial_merge!(op::Union{typeof(+), typeof(-), MA.AddSubMul}, p::Polynomial{T,TT}, args...) where {T,TT}
    get1 = let p=p
        i -> p.terms[i]
    end
    set = let p=p
        (i, t) -> begin
            p.terms[i] = t
        end
    end
    push = let p=p
        t -> push!(p.terms, t)
    end
    resize = let p=p
        (n) -> resize!(p.terms, n)
    end
    # We can modify the coefficient since it's the result of `combine`.
    keep = let p=p
        i -> begin
            if i isa Int
                !MA.iszero!!(coefficient(p.terms[i]))
            else
                !MA.iszero!!(coefficient(i))
            end
        end
    end
    __polynomial_merge!(op, p, get1, set, push, resize, keep, args...)
    return p
end

function MA.operate!(op::Union{typeof(+), typeof(-)}, p::Polynomial, q::Union{AbstractTermLike,Polynomial})
    return _polynomial_merge!(op, p, q)
end

function MA.operate!(op::MA.AddSubMul, p::Polynomial, q::Polynomial, args::AbstractTermLike...)
    return _polynomial_merge!(op, p, q, args...)
end

function MA.operate!(op::MA.AddSubMul, p::Polynomial, t::AbstractTermLike, q::Polynomial)
    return _polynomial_merge!(op, p, t, q)
end

function MA.buffer_for(op::MA.AddSubMul, ::Type{<:Polynomial{S}}, ::Type{<:AbstractTermLike{T}}, ::Type{<:Polynomial{U}}) where {S,T,U}
    return MA.buffer_for(op, S, T, U)
end
function MA.buffered_operate!(buffer, op::MA.AddSubMul, p::Polynomial, t::AbstractTermLike, q::Polynomial)
    return _polynomial_merge!(op, p, t, q, buffer)
end

function MA.operate_to!(output::Polynomial, ::typeof(*), p::Polynomial, q::Polynomial)
    empty!(output.terms)
    mul_to_terms!(output.terms, p, q)
    sort!(output.terms, lt=(<))
    uniqterms!(output.terms)
    return output
end
function MA.operate!(::typeof(*), p::Polynomial, q::Polynomial)
    if iszero(q)
        return MA.operate!(zero, p)
    elseif nterms(q) == 1
        return MA.operate!(*, p, leading_term(q))
    else
        return MA.operate_to!(p, *, MA.mutable_copy(p), q)
    end
end
function MA.operate!(::typeof(*), p::Polynomial, t::AbstractTermLike)
    for i in eachindex(p.terms)
        p.terms[i] = MA.operate!!(*, p.terms[i], t)
    end
    return p
end
function map_exponents!(f, p::Polynomial, m::AbstractMonomialLike)
    for i in eachindex(p.terms)
        t = p.terms[i]
        p.terms[i] = term(coefficient(t), map_exponents(f, monomial(t), m))
    end
    return p
end
function map_exponents(f, p::Polynomial, m::AbstractMonomialLike)
    P = MA.promote_operation(*, typeof(p), typeof(m))
    return map_exponents!(f, MA.mutable_copy(convert(P, p)), m)
end

function MA.operate!(::typeof(zero), p::Polynomial)
    empty!(p.terms)
    return p
end
function MA.operate!(::typeof(one), p::Polynomial{T}) where T
    if isempty(p.terms)
        push!(p.terms, constant_term(one(T), p))
    else
        t = p.terms[1]
        p.terms[1] = Term(MA.one!!(coefficient(t)), constant_monomial(t))
        resize!(p.terms, 1)
    end
    return p
end

function MA.operate!(::typeof(remove_leading_term), p::Polynomial)
    pop!(p.terms)
    return p
end

function MA.operate!(::typeof(unsafe_restore_leading_term), p::Polynomial, t::AbstractTermLike)
    # We don't need to copy the coefficient of `t`, this is why this function is called `unsafe`
    push!(p.terms, t)
    return p
end

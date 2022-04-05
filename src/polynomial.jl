export Polynomial

struct Polynomial{CoeffType,T<:AbstractTerm{CoeffType},V<:AbstractVector{T}} <: AbstractPolynomial{CoeffType}
    terms::V

    Polynomial{C, T, V}(terms::AbstractVector{T}) where {C, T, V} = new{C, T, V}(terms)
end

const VectorPolynomial{C,T} = Polynomial{C,T,Vector{T}}

termtype(::Type{<:Polynomial{C,T}}) where {C,T} = T
terms(p::Polynomial) = p.terms
constantmonomial(::Union{Polynomial{C,TT},Type{Polynomial{C,TT}}}) where {C,TT} = constantmonomial(TT)
function convertconstant(::Type{Polynomial{C,TT,Vector{TT}}}, α) where {C,TT}
    t = term(convert(C, α), constantmonomial(TT))
    return Polynomial{C,TT,Vector{TT}}([t])
end
Base.convert(PT::Type{Polynomial{C,TT,V}}, p::Polynomial) where {C,TT,V} = PT(convert(V, p.terms))

_change_eltype(::Type{<:Vector}, ::Type{T}) where T = Vector{T}
function polynomialtype(::Type{Polynomial{C, T, V}}, ::Type{NewC}) where {C, T, V, NewC}
    NewT = termtype(T, NewC)
    Polynomial{NewC, NewT, _change_eltype(V, NewT)}
end

(p::Polynomial)(s...) = substitute(Eval(), p, s)

function Base.one(::Type{P}) where {C,TT,P<:Polynomial{C,TT}}
    return convert(P, one(TT))
end
Base.one(p::Polynomial) = one(typeof(p))

Base.zero(::Type{Polynomial{C,T,A}}) where {C,T,A} = Polynomial{C,T,A}(A())
Base.zero(t::Polynomial) = zero(typeof(t))

jointerms(terms1::AbstractArray{<:Term}, terms2::AbstractArray{<:Term}) = Sequences.mergesorted(terms1, terms2, compare, combine)
function jointerms!(output::AbstractArray{<:Term}, terms1::AbstractArray{<:Term}, terms2::AbstractArray{<:Term})
    resize!(output, length(terms1) + length(terms2))
    Sequences.mergesorted!(output, terms1, terms2, compare, combine)
end

Base.:(+)(p1::Polynomial, p2::Polynomial) = polynomial!(jointerms(terms(p1), terms(p2)), SortedUniqState())
Base.:(+)(p1::Polynomial, p2::Polynomial{<:LinearAlgebra.UniformScaling}) = p1 + mapcoefficientsnz(J -> J.λ, p2)
Base.:(+)(p1::Polynomial{<:LinearAlgebra.UniformScaling}, p2::Polynomial) = mapcoefficientsnz(J -> J.λ, p1) + p2
Base.:(+)(p1::Polynomial{<:LinearAlgebra.UniformScaling}, p2::Polynomial{<:LinearAlgebra.UniformScaling}) = mapcoefficientsnz(J -> J.λ, p1) + p2
function MA.operate_to!(result::Polynomial, ::typeof(+), p1::Polynomial, p2::Polynomial)
    if result === p1 || result === p2
        error("Cannot call `operate_to!(output, +, p, q)` with `output` equal to `p` or `q`, call `operate!` instead.")
    end
    jointerms!(result.terms, terms(p1), terms(p2))
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
Base.:(-)(p1::Polynomial, p2::Polynomial) = polynomial!(jointerms(terms(p1), (-).(terms(p2))))
Base.:(-)(p1::Polynomial, p2::Polynomial{<:LinearAlgebra.UniformScaling}) = p1 - mapcoefficientsnz(J -> J.λ, p2)
Base.:(-)(p1::Polynomial{<:LinearAlgebra.UniformScaling}, p2::Polynomial) = mapcoefficientsnz(J -> J.λ, p1) - p2
Base.:(-)(p1::Polynomial{<:LinearAlgebra.UniformScaling}, p2::Polynomial{<:LinearAlgebra.UniformScaling}) = mapcoefficientsnz(J -> J.λ, p1) - p2

LinearAlgebra.adjoint(x::Polynomial) = polynomial!(adjoint.(terms(x)))

function mapcoefficients(f::Function, p::Polynomial; nonzero = false)
    terms = map(p.terms) do term
        mapcoefficients(f, term)
    end
    if !nonzero
        filter!(!iszero, terms)
    end
    return polynomial!(terms)
end
function mapcoefficients!(f::Function, p::Polynomial; nonzero = false)
    for i in eachindex(p.terms)
        t = p.terms[i]
        p.terms[i] = Term(f(coefficient(t)), monomial(t))
    end
    if !nonzero
        filter!(!iszero, p.terms)
    end
    return p
end

function mapcoefficients_to!(output::Polynomial, f::Function, p::Polynomial; nonzero = false)
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
function mapcoefficients_to!(output::Polynomial, f::Function, p::AbstractPolynomialLike; nonzero = false)
    return mapcoefficients_to!(output, f, polynomial(p); nonzero = false)
end

# The polynomials can be mutated.
MA.mutability(::Type{<:VectorPolynomial}) = MA.IsMutable()

# Terms are not mutable under the MutableArithmetics API
function MA.mutable_copy(p::VectorPolynomial{C,TT}) where {C,TT}
    return VectorPolynomial{C,TT}([TT(MA.copy_if_mutable(coefficien(term)), monomial(term)) for term in terms(p)])
end
Base.copy(p::VectorPolynomial) = MA.mutable_copy(p)

function grlex end
function MA.operate!(op::Union{typeof(+), typeof(-)}, p::Polynomial{T,TT}, q::Polynomial) where {T,TT}
    get1(i) = p.terms[i]
    function get2(i)
        t = q.terms[i]
        TT(MA.scaling_convert(T, MA.operate(op, coefficient(t))), monomial(t))
    end
    set(i, t::TT) = p.terms[i] = t
    push(t::TT) = push!(p.terms, t)
    compare_monomials(t::TT, j::Int) = grlex(monomial(q.terms[j]), monomial(t))
    compare_monomials(i::Int, j::Int) = compare_monomials(get1(i), j)
    combine(i::Int, j::Int) = p.terms[i] = Term(MA.operate!!(op, coefficient(p.terms[i]), coefficient(q.terms[j])), monomial(p.terms[i]))
    combine(t::TT, j::Int) = TT(MA.operate!!(op, coefficient(t), coefficient(q.terms[j])), monomial(t))
    resize(n) = resize!(p.terms, n)
    # We can modify the coefficient since it's the result of `combine`.
    keep(t::Term) = !MA.iszero!!(coefficient(t))
    keep(i::Int) = !MA.iszero!!(coefficient(p.terms[i]))
    polynomial_merge!(
        nterms(p), nterms(q), get1, get2, set, push,
        compare_monomials, combine, keep, resize
    )
    return p
end
function MA.operate_to!(output::Polynomial, ::typeof(*), p::Polynomial, q::Polynomial)
    empty!(output.terms)
    mul_to_terms!(output.terms, p, q)
    sort!(output.terms, lt=(>))
    uniqterms!(output.terms)
    return output
end
function MA.operate!(::typeof(*), p::Polynomial, q::Polynomial)
    return MA.operate_to!(p, *, MA.mutable_copy(p), q)
end

function MA.operate!(::typeof(zero), p::Polynomial)
    empty!(p.terms)
    return p
end
function MA.operate!(::typeof(one), p::Polynomial{T}) where T
    if isempty(p.terms)
        push!(p.terms, constantterm(one(T), p))
    else
        t = p.terms[1]
        p.terms[1] = Term(MA.one!!(coefficient(t)), constantmonomial(t))
        resize!(p.terms, 1)
    end
    return p
end

function MA.operate!(::typeof(removeleadingterm), p::Polynomial)
    deleteat!(p.terms, 1)
    return p
end

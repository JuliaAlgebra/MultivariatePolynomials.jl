export NCTerm, NCVecPolynomial, NCMatPolynomial, NCSOSDecomposition, getmat, monomials, removemonomials, NCTermType

abstract NCTermType{T} <: NCPolyType
zero(::Type{NCPolyType}) = zero(NCVecPolynomial{Int})
one(::Type{NCPolyType}) = one(NCVecPolynomial{Int})
zero(p::NCPolyType) = zero(NCPolyType)
one(p::NCPolyType) = one(NCPolyType)

zero{T}(t::NCTermType{T}) = NCVecPolynomial(T[], NCMonomialVector(vars(t), Vector{Vector{Int}}()))
zero{T<:NCTermType}(::Type{T}) = NCVecPolynomial(eltype(T)[], NCMonomialVector())
one{T}(t::NCTermType{T}) = NCVecPolynomial([one(T)], NCMonomialVector(vars(t), [zeros(Int, length(vars(t)))]))
one{T<:NCTermType}(::Type{T}) = NCVecPolynomial([one(eltype(T))], NCMonomialVector(PolyVar[], [Int[]]))

abstract NCTermContainer{T} <: NCTermType{T}
eltype{T}(::Type{NCTermContainer{T}}) = T
eltype{T}(::Type{NCTermType{T}}) = T

type NCTerm{T} <: NCTermContainer{T}
    α::T
    x::NCMonomial
end
NCTerm(t::NCTerm) = t
NCTerm(x::NCMonomial) = NCTerm{Int}(x)
NCTerm(x::NCPolyVar) = NCTerm(NCMonomial(x))
NCTermContainer{T<:Union{NCMonomial,NCPolyVar}}(x::T) = NCTerm(x)
NCTermContainer{T}(t::NCTermContainer{T}) = t
NCTerm{T}(α::T) = NCTerm{T}(α, NCMonomial())
Base.convert{T}(::Type{NCTerm{T}}, t::NCTerm{T}) = t
Base.convert{T}(::Type{NCTerm{T}}, t::NCTerm) = NCTerm{T}(T(t.α), t.x)
Base.convert{T}(::Type{NCTerm{T}}, x::NCMonomial) = NCTerm{T}(one(T), x)
Base.convert{T}(::Type{NCTerm{T}}, x::NCPolyVar) = NCTerm{T}(NCMonomial(x))
Base.convert{T}(::Type{NCTerm{T}}, α) = NCTerm(T(α))

Base.convert{T}(::Type{NCTermContainer{T}}, x::Union{NCMonomial,NCPolyVar}) = NCTerm{T}(x)
Base.convert{T}(::Type{NCTermContainer{T}}, α::T) = NCTerm{T}(α, NCMonomial())
Base.convert{S,T}(::Type{NCTermContainer{T}}, α::S) = NCTermContainer{T}(T(α))
NCTermContainer{T}(α::T) = NCTermContainer{T}(α)

Base.convert{T}(::Type{NCTermContainer{T}}, t::NCTerm) = NCTerm{T}(t)
Base.convert{T<:NCTermContainer}(::Type{T}, t::NCTerm) = NCTerm{eltype(T)}(t)

zero{T}(t::NCTerm{T}) = NCTerm(zero(T), t.x)
zero{T}(::Type{NCTerm{T}}) = NCTerm(zero(T), Monomial())
one{T}(t::NCTerm{T}) = NCTerm(one(T), Monomial(t.x.vars, zeros(Int, length(t.x.vars))))
one{T}(::Type{NCTerm{T}}) = NCTerm(one(T), Monomial())

# Invariant:
# a and x might be empty: meaning it is the zero polynomial
# a does not contain any zeros
# x is increasing in the monomial order (i.e. grlex)
type NCVecPolynomial{T} <: NCTermContainer{T}
    a::Vector{T}
    x::NCMonomialVector

    function NCVecPolynomial(a::Vector{T}, x::NCMonomialVector)
        if length(a) != length(x)
            error("There should be as many coefficient than monomials")
        end
        zeroidx = Int[]
        for (i,α) in enumerate(a)
            if iszero(α)
                push!(zeroidx, i)
            end
        end
        if !isempty(zeroidx)
            isnz = ones(Bool, length(a))
            isnz[zeroidx] = false
            nzidx = find(isnz)
            a = a[nzidx]
            x = x[nzidx]
        end
        new(a, x)
    end
end
(::Type{NCVecPolynomial{T}}){T}(a::Vector{T}, x::Vector) = NCVecPolynomial{T}(a, NCMonomialVector(x))
(::Type{NCVecPolynomial{T}}){S,T}(a::Vector{S}, x::Vector) = NCVecPolynomial{T}(Vector{T}(a), NCMonomialVector(x))

Base.copy{T}(p::NCVecPolynomial{T}) = NCVecPolynomial{T}(copy(p.a), copy(p.x))

NCVecPolynomial{T}(a::Vector{T}, x::NCMonomialVector) = NCVecPolynomial{T}(a, x)
NCVecPolynomial(a::Vector, x::Vector) = NCVecPolynomial(a, NCMonomialVector(x))

NCVecPolynomial(x) = NCVecPolynomial(NCTerm(x))
NCVecPolynomial{T}(p::NCVecPolynomial{T}) = p
NCVecPolynomial{T}(t::NCTerm{T}) = NCVecPolynomial{T}([t.α], [t.x])
Base.convert{T}(::Type{NCVecPolynomial{T}}, x) = NCVecPolynomial(NCTerm{T}(x))
Base.convert{T}(::Type{NCVecPolynomial{T}}, t::NCTerm) = NCVecPolynomial{T}([T(t.α)], [t.x])
Base.convert{T}(::Type{NCVecPolynomial{T}}, p::NCVecPolynomial{T}) = p
Base.convert{S,T}(::Type{NCVecPolynomial{T}}, p::NCVecPolynomial{S}) = NCVecPolynomial(Vector{T}(p.a), p.x)

Base.convert{T}(::Type{NCTermContainer{T}}, p::NCVecPolynomial) = NCVecPolynomial{T}(p)

function (::Type{NCVecPolynomial{T}}){T}(f::Function, x::NCMonomialVector)
    a = T[f(i) for i in 1:length(x)]
    NCVecPolynomial{T}(a, x)
end
(::Type{NCVecPolynomial{T}}){T}(f::Function, x::Vector) = NCVecPolynomial{T}(f, NCMonomialVector(x))

function vecpolynomialclean{T}(vars::Vector{NCPolyVar}, adup::Vector{T}, Zdup::Vector{Vector{Int}})
    a, Z = removedups(adup, Zdup)
    NCVecPolynomial(a, NCMonomialVector(vars, Z))
end

eltype{T}(::Type{NCVecPolynomial{T}}) = T
getindex(p::NCVecPolynomial, I::Int) = NCTerm(p.a[I[1]], p.x[I[1]])
removemonomials(p::NCVecPolynomial, x::Vector) = removemonomials(p, NCMonomialVector(x))

type NCMatPolynomial{T} <: NCTermType{T}
    Q::Vector{T}
    x::NCMonomialVector
end

function (::Type{NCMatPolynomial{T}}){T}(f::Function, x::NCMonomialVector)
    NCMatPolynomial{T}(trimat(f, length(x)), x)
end
(::Type{NCMatPolynomial{T}}){T}(f::Function, x::Vector) = NCMatPolynomial{T}(f, NCMonomialVector(x))

function NCMatPolynomial{T}(Q::Matrix{T}, x::NCMonomialVector)
    NCMatPolynomial{T}((i,j) -> Q[i,j], x)
end
NCMatPolynomial(Q::Matrix, x::Vector) = NCMatPolynomial(Q, NCMonomialVector(x))

function NCVecPolynomial{T}(p::NCMatPolynomial{T})
    if isempty(p.Q)
        zero(NCVecPolynomial{T})
    else
        n = length(p.x)
        N = trimap(n, n, n)
        Z = Vector{Vector{Int}}(N)
        U = typeof(2*p.Q[1] + p.Q[1])
        a = Vector{U}(N)
        for i in 1:n
            for j in i:n
                k = trimap(i, j, n)
                Z[k] = p.x.Z[i] + p.x.Z[j] # FIXME this is wrong for NC
                if i == j
                    a[k] = p.Q[k]
                else
                    a[k] = 2*p.Q[k]
                end
            end
        end
        vecpolynomialclean(p.x.vars, a, Z)
    end
end
NCTermContainer(p::NCMatPolynomial) = NCVecPolynomial(p)

type NCSOSDecomposition{T} <: NCTermType{T}
    ps::Vector{NCVecPolynomial{T}}
    function NCSOSDecomposition(ps::Vector{NCVecPolynomial{T}})
        new(ps)
    end
end
function (::Type{NCSOSDecomposition{T}}){T}(ps::Vector)
    NCSOSDecomposition(Vector{NCVecPolynomial{T}}(ps))
end
function NCSOSDecomposition(ps::Vector)
    T = reduce(promote_type, Int, map(eltype, ps))
    NCSOSDecomposition{T}(ps)
end

function NCSOSDecomposition{T}(p::NCMatPolynomial{T})
    n = length(p.x)
    # TODO LDL^T factorization for SDP is missing in Julia
    # it would be nice to have though
    A = getmat(p)
    Q = chol(A)
    ps = [NCVecPolynomial(Q[i,:], p.x) for i in 1:n]
    NCSOSDecomposition(ps)
end

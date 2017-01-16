export Term, VecPolynomial, MatPolynomial, SOSDecomposition, TermType

abstract TermType{T} <: PolyType
zero(::Type{PolyType}) = zero(VecPolynomial{Int})
one(::Type{PolyType}) = one(VecPolynomial{Int})
zero(p::PolyType) = zero(PolyType)
one(p::PolyType) = one(PolyType)

zero{T}(t::TermType{T}) = VecPolynomial(T[], MonomialVector(vars(t), Vector{Vector{Int}}()))
zero{T<:TermType}(::Type{T}) = VecPolynomial(eltype(T)[], MonomialVector())
one{T}(t::TermType{T}) = VecPolynomial([one(T)], MonomialVector(vars(t), [zeros(Int, length(vars(t)))]))
one{T<:TermType}(::Type{T}) = VecPolynomial([one(eltype(T))], MonomialVector(PolyVar[], [Int[]]))

abstract TermContainer{T} <: TermType{T}
eltype{T}(::Type{TermContainer{T}}) = T
eltype{T}(::Type{TermType{T}}) = T

type Term{T} <: TermContainer{T}
    α::T
    x::Monomial
end
Term(t::Term) = t
Term(x::Monomial) = Term{Int}(x)
Term(x::PolyVar) = Term(Monomial(x))
TermContainer{T<:Union{Monomial,PolyVar}}(x::T) = Term(x)
TermContainer{T}(t::TermContainer{T}) = t
Term{T}(α::T) = Term{T}(α, Monomial())
Base.convert{T}(::Type{Term{T}}, t::Term{T}) = t
Base.convert{T}(::Type{Term{T}}, t::Term) = Term{T}(T(t.α), t.x)
Base.convert{T}(::Type{Term{T}}, x::Monomial) = Term{T}(one(T), x)
Base.convert{T}(::Type{Term{T}}, x::PolyVar) = Term{T}(Monomial(x))
Base.convert{T}(::Type{Term{T}}, α) = Term(T(α))

Base.convert{T}(::Type{TermContainer{T}}, x::Union{Monomial,PolyVar}) = Term{T}(x)
Base.convert{T}(::Type{TermContainer{T}}, α::T) = Term{T}(α, Monomial())
Base.convert{S,T}(::Type{TermContainer{T}}, α::S) = TermContainer{T}(T(α))
TermContainer{T}(α::T) = TermContainer{T}(α)

Base.convert{T}(::Type{TermContainer{T}}, t::Term) = Term{T}(t)
Base.convert{T<:TermContainer}(::Type{T}, t::Term) = Term{eltype(T)}(t)

zero{T}(t::Term{T}) = Term(zero(T), t.x)
zero{T}(::Type{Term{T}}) = Term(zero(T), Monomial())
one{T}(t::Term{T}) = Term(one(T), Monomial(t.x.vars, zeros(Int, length(t.x.vars))))
one{T}(::Type{Term{T}}) = Term(one(T), Monomial())

# Invariant:
# a and x might be empty: meaning it is the zero polynomial
# a does not contain any zeros
# x is increasing in the monomial order (i.e. grlex)
type VecPolynomial{T} <: TermContainer{T}
    a::Vector{T}
    x::MonomialVector

    function VecPolynomial(a::Vector{T}, x::MonomialVector)
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
(::Type{VecPolynomial{T}}){T}(a::Vector{T}, x::Vector) = VecPolynomial{T}(a, MonomialVector(x))
(::Type{VecPolynomial{T}}){S,T}(a::Vector{S}, x::Vector) = VecPolynomial{T}(Vector{T}(a), MonomialVector(x))

Base.copy{T}(p::VecPolynomial{T}) = VecPolynomial{T}(copy(p.a), copy(p.x))

VecPolynomial{T}(a::Vector{T}, x::MonomialVector) = VecPolynomial{T}(a, x)
VecPolynomial(a::Vector, x::Vector) = VecPolynomial(a, MonomialVector(x))

VecPolynomial(x) = VecPolynomial(Term(x))
VecPolynomial{T}(p::VecPolynomial{T}) = p
VecPolynomial{T}(t::Term{T}) = VecPolynomial{T}([t.α], [t.x])
Base.convert{T}(::Type{VecPolynomial{T}}, x) = VecPolynomial(Term{T}(x))
Base.convert{T}(::Type{VecPolynomial{T}}, t::Term) = VecPolynomial{T}([T(t.α)], [t.x])
Base.convert{T}(::Type{VecPolynomial{T}}, p::VecPolynomial{T}) = p
Base.convert{S,T}(::Type{VecPolynomial{T}}, p::VecPolynomial{S}) = VecPolynomial(Vector{T}(p.a), p.x)

Base.convert{T}(::Type{TermContainer{T}}, p::VecPolynomial) = VecPolynomial{T}(p)

function (::Type{VecPolynomial{T}}){T}(f::Function, x::MonomialVector)
    a = T[f(i) for i in 1:length(x)]
    VecPolynomial{T}(a, x)
end
(::Type{VecPolynomial{T}}){T}(f::Function, x::Vector) = VecPolynomial{T}(f, MonomialVector(x))

function vecpolynomialclean{T}(vars::Vector{PolyVar}, adup::Vector{T}, Zdup::Vector{Vector{Int}})
    a, Z = removedups(adup, Zdup)
    VecPolynomial(a, MonomialVector(vars, Z))
end

eltype{T}(::Type{VecPolynomial{T}}) = T
getindex(p::VecPolynomial, I::Int) = Term(p.a[I[1]], p.x[I[1]])
removemonomials(p::VecPolynomial, x::Vector) = removemonomials(p, MonomialVector(x))

type MatPolynomial{T} <: TermType{T}
    Q::Vector{T}
    x::MonomialVector
end

function (::Type{MatPolynomial{T}}){T}(f::Function, x::MonomialVector)
    MatPolynomial{T}(trimat(f, length(x)), x)
end
(::Type{MatPolynomial{T}}){T}(f::Function, x::Vector) = MatPolynomial{T}(f, MonomialVector(x))

function MatPolynomial{T}(Q::Matrix{T}, x::MonomialVector)
    MatPolynomial{T}((i,j) -> Q[i,j], x)
end
MatPolynomial(Q::Matrix, x::Vector) = MatPolynomial(Q, MonomialVector(x))

function VecPolynomial{T}(p::MatPolynomial{T})
    if isempty(p.Q)
        zero(VecPolynomial{T})
    else
        n = length(p.x)
        N = trimap(n, n, n)
        Z = Vector{Vector{Int}}(N)
        U = typeof(2*p.Q[1] + p.Q[1])
        a = Vector{U}(N)
        for i in 1:n
            for j in i:n
                k = trimap(i, j, n)
                Z[k] = p.x.Z[i] + p.x.Z[j]
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
TermContainer(p::MatPolynomial) = VecPolynomial(p)

type SOSDecomposition{T} <: TermType{T}
    ps::Vector{VecPolynomial{T}}
    function SOSDecomposition(ps::Vector{VecPolynomial{T}})
        new(ps)
    end
end
function (::Type{SOSDecomposition{T}}){T}(ps::Vector)
    SOSDecomposition(Vector{VecPolynomial{T}}(ps))
end
function SOSDecomposition(ps::Vector)
    T = reduce(promote_type, Int, map(eltype, ps))
    SOSDecomposition{T}(ps)
end

function SOSDecomposition{T}(p::MatPolynomial{T})
    n = length(p.x)
    # TODO LDL^T factorization for SDP is missing in Julia
    # it would be nice to have though
    A = getmat(p)
    Q = chol(A)
    ps = [VecPolynomial(Q[i,:], p.x) for i in 1:n]
    SOSDecomposition(ps)
end

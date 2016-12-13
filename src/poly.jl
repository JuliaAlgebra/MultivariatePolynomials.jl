export Term, VecPolynomial, MatPolynomial, SOSDecomposition, getmat, monomials, removemonomials, TermType
import Base.eltype, Base.zero, Base.one

abstract TermType{T} <: PolyType
zero{T}(t::TermType{T}) = VecPolynomial(T[], MonomialVector(vars(t), Vector{Vector{Int}}()))
zero(::Type{PolyType}) = zero(VecPolynomial{Int})
zero{T<:TermType}(::Type{T}) = VecPolynomial(eltype(T)[], MonomialVector(PolyVar[], Vector{Vector{Int}}()))
one{T}(t::TermType{T}) = Term(one(T), Monomial(vars(t), zeros(Int, length(vars(t)))))
one{T<:TermType}(::Type{T}) = Term(one(eltype(T)), Monomial(PolyVar[], Int[]))

abstract TermContainer{T} <: TermType{T}

eltype{T<:TermType}(::Type{T}) = T.parameters[1]
eltype{T}(p::TermType{T}) = T

# Invariant:
# α is nonzero (otherwise, just keep zero(T) and drop the monomial x)
type Term{T} <: TermContainer{T}
    α::T
    x::Monomial
end
Term(t::Term) = t
Term(x::Monomial) = Term{Int}(x)
Term(x::PolyVar) = Term(Monomial(x))
TermContainer{T<:Union{Monomial,PolyVar}}(x::T) = Term(x)
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

Base.convert{T<:TermContainer}(::Type{T}, t::Term) = Term{eltype(T)}(t)
Base.convert(::Type{Any}, t::Term) = t

vars(t::Term) = vars(t.x)

length(::Term) = 1
isempty(::Term) = false
start(::Term) = false
done(::Term, state) = state
next(t::Term, state) = (t, true)
getindex(t::Term, I::Int) = t

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
Base.convert(::Type{Any}, p::VecPolynomial) = p

Base.convert(::Type{PolyType}, p::TermContainer) = p
function Base.convert{S}(::Type{S}, p::TermContainer)
    s = zero(S)
    for t in p
        if sum(abs(t.x.z)) > 0
            # The polynomial is not constant
            throw(InexactError())
        end
        s += S(t.α)
    end
    s
end

function (::Type{VecPolynomial{T}}){T}(f::Function, x::MonomialVector)
    a = T[f(i) for i in 1:length(x)]
    VecPolynomial{T}(a, x)
end
(::Type{VecPolynomial{T}}){T}(f::Function, x::Vector) = VecPolynomial{T}(f, MonomialVector(x))

function vecpolynomialclean{T}(vars::Vector{PolyVar}, adup::Vector{T}, Zdup::Vector{Vector{Int}})
    σ = sortperm(Zdup, rev=true)
    Z = Vector{Vector{Int}}()
    a = Vector{T}()
    i = 0
    j = 1
    while j <= length(adup)
        k = σ[j]
        if j == 1 || Zdup[k] != Zdup[σ[j-1]]
            push!(Z, Zdup[k])
            push!(a, adup[k])
            i += 1
        else
            a[i] += adup[k]
        end
        j += 1
    end
    VecPolynomial(a, MonomialVector(vars, Z))
end

vars(p::VecPolynomial) = vars(p.x)

length(p::VecPolynomial) = length(p.a)
isempty(p::VecPolynomial) = length(p) > 0
start(::VecPolynomial) = 1
done(p::VecPolynomial, state) = length(p) < state
next(p::VecPolynomial, state) = (p[state], state+1)
getindex(p::VecPolynomial, I::Int) = Term(p.a[I[1]], p.x[I[1]])

monomials(p::VecPolynomial) = p.x
function removemonomials(p::VecPolynomial, x::MonomialVector)
    # use the fact that monomials are sorted to do this O(n) instead of O(n^2)
    j = 1
    I = Int[]
    for (i,t) in enumerate(p)
        while j <= length(x) && x[j] < t.x
            j += 1
        end
        if j > length(x) || x[j] != t.x
            push!(I, i)
        end
    end
    VecPolynomial(p.a[I], p.x[I])
end
removemonomials(p::VecPolynomial, x::Vector) = removemonomials(p, MonomialVector(x))

type MatPolynomial{T} <: TermType{T}
    Q::Vector{T}
    x::MonomialVector
end

function trimap(i, j, n)
    div(n*(n+1), 2) - div((n-i+1)*(n-i+2), 2) + j-i+1
end

function (::Type{MatPolynomial{T}}){T}(f::Function, x::MonomialVector)
    n = length(x)
    Q = Vector{T}(trimap(n, n, n))
    for i in 1:n
        for j in i:n
            Q[trimap(i,j,n)] = f(i,j)
        end
    end
    MatPolynomial{T}(Q, x)
end
(::Type{MatPolynomial{T}}){T}(f::Function, x::Vector) = MatPolynomial{T}(f, MonomialVector(x))

function MatPolynomial{T}(Q::Matrix{T}, x::MonomialVector)
    MatPolynomial{T}((i,j) -> Q[i,j], x)
end
MatPolynomial(Q::Matrix, x::Vector) = MatPolynomial(Q, MonomialVector(x))

function getindex(p::MatPolynomial, I::NTuple{2,Int})
    i, j = I
    if i < j
        i, j = (j, i)
    end
    n = length(p.x)
    p.Q[trimap(i,j,n)]
end

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
length(p::SOSDecomposition) = length(p.ps)
isempty(p::SOSDecomposition) = isempty(p.ps)
start(p::SOSDecomposition) = start(p.ps)
done(p::SOSDecomposition, state) = done(p.ps, state)
next(p::SOSDecomposition, state) = next(p.ps, state)

function getmat{T}(p::MatPolynomial{T})
    n = length(p.x)
    A = Matrix{T}(n, n)
    for i in 1:n, j in i:n
        A[j,i] = A[i,j] = p.Q[trimap(i,j,n)]
    end
    A
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

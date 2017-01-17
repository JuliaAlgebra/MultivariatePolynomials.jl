export Term, VecPolynomial, MatPolynomial, SOSDecomposition, TermType, getmat, monomials, removemonomials

abstract TermType{C, T} <: PolyType{C}
eltype{C, T}(::Type{TermType{C, T}}) = T
eltype{C, T}(p::TermType{C, T}) = T
zero{C, T}(t::TermType{C, T}) = VecPolynomial(T[], MonomialVector{C}(vars(t), Vector{Vector{Int}}()))
#zero{T<:TermType}(::Type{T}) = VecPolynomial(eltype(T)[], MonomialVector{iscomm(T)}())
one{C, T}(t::TermType{C, T}) = VecPolynomial([one(T)], MonomialVector{C}(vars(t), [zeros(Int, length(vars(t)))]))
#one{T<:TermType}(::Type{T}) = VecPolynomial([one(eltype(T))], MonomialVector{iscomm(T)}(PolyVar[], [Int[]]))

abstract TermContainer{C, T} <: TermType{C, T}
eltype{C, T}(::Type{TermContainer{C, T}}) = T
zero{C, T}(::Type{TermContainer{C, T}}) = zero(VecPolynomial{C, T})
one{C, T}(::Type{TermContainer{C, T}}) = one(VecPolynomial{C, T})

type Term{C, T} <: TermContainer{C, T}
    α::T
    x::Monomial{C}
end
iscomm{C, T}(::Type{Term{C, T}}) = C
Term(t::Term) = t
(::Type{Term{C}}){C}(x::Monomial{C}) = Term{C, Int}(x)
(::Type{Term{C}}){C}(x::PolyVar{C}) = Term{C}(Monomial{C}(x))
Term{C}(x::Monomial{C}) = Term{C}(x)
Term{C}(x::PolyVar{C}) = Term{C}(x)
(::Type{TermContainer{C}}){C}(x::PolyVar{C}) = Term(x)
(::Type{TermContainer{C}}){C}(x::Monomial{C}) = Term(x)
(::Type{TermContainer{C}}){C}(t::TermContainer{C}) = t
TermContainer(x::PolyVar) = Term(x)
TermContainer(x::Monomial) = Term(x)
TermContainer(t::TermContainer) = t
(::Type{Term{C}}){C, T}(α::T) = Term{C, T}(α, Monomial{C}())
Base.convert{C, T}(::Type{Term{C, T}}, t::Term{C, T}) = t
Base.convert{C, T}(::Type{Term{C, T}}, t::Term{C}) = Term{C, T}(T(t.α), t.x)
Base.convert{C, T}(::Type{Term{C, T}}, x::Monomial{C}) = Term{C, T}(one(T), x)
Base.convert{C, T}(::Type{Term{C, T}}, x::PolyVar{C}) = Term{C, T}(Monomial{C}(x))
Base.convert{C, T}(::Type{Term{C, T}}, α) = Term{C}(T(α))

Base.convert{C, T}(::Type{TermContainer{C, T}}, x::Union{Monomial{C},PolyVar{C}}) = Term{C, T}(x)
Base.convert{C, T}(::Type{TermContainer{C, T}}, α::T) = Term{C, T}(α, Monomial{C}())
Base.convert{C, S, T}(::Type{TermContainer{C, T}}, α::S) = TermContainer{C, T}(T(α))
(::Type{TermContainer{C}}){C, T}(α::T) = TermContainer{C, T}(α)

Base.convert{C, T}(::Type{TermContainer{C, T}}, t::Term{C}) = Term{C, T}(t)
Base.convert{T<:TermContainer}(::Type{T}, t::Term) = Term{iscomm(T), eltype(T)}(t)

zero{C, T}(t::Term{C, T}) = Term{C, T}(zero(T), t.x)
zero{C, T}(::Type{Term{C, T}}) = Term{C, T}(zero(T), Monomial{C}())
one{C, T}(t::Term{C, T}) = Term{C, T}(one(T), Monomial{C}(t.x.vars, zeros(Int, length(t.x.vars))))
one{C, T}(::Type{Term{C, T}}) = Term{C, T}(one(T), Monomial{C}())

Base.convert(::Type{Any}, t::Term) = t
Base.copy{T<:Term}(t::T) = T(copy(t.α), copy(t.x))

vars(t::Term) = vars(t.x)

eltype{T}(::Type{Term{T}}) = T
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
type VecPolynomial{C, T} <: TermContainer{C, T}
    a::Vector{T}
    x::MonomialVector{C}

    function VecPolynomial(a::Vector{T}, x::MonomialVector{C})
        if length(a) != length(x)
            throw(ArgumentError("There should be as many coefficient than monomials"))
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
iscomm{C, T}(::Type{VecPolynomial{C, T}}) = C
(::Type{VecPolynomial{C, T}}){C, T}(a::Vector{T}, x::Vector) = VecPolynomial{C, T}(a, MonomialVector{C}(x))
(::Type{VecPolynomial{C, T}}){C, S, T}(a::Vector{S}, x::Vector) = VecPolynomial{C, T}(Vector{T}(a), MonomialVector{C}(x))

Base.copy{C, T}(p::VecPolynomial{C, T}) = VecPolynomial{C, T}(copy(p.a), copy(p.x))
zero{C, T}(::Type{VecPolynomial{C, T}}) = VecPolynomial(T[], MonomialVector{C}())
one{C, T}(::Type{VecPolynomial{C, T}}) = VecPolynomial([one(T)], MonomialVector{C}(PolyVar{C}[], [Int[]]))

VecPolynomial{C, T}(a::Vector{T}, x::MonomialVector{C}) = VecPolynomial{C, T}(a, x)
(::Type{VecPolynomial{C}}){C, T}(a::Vector{T}, x::MonomialVector{C}) = VecPolynomial{C, T}(a, x)
(::Type{VecPolynomial{C}}){C}(a::Vector, x::Vector) = VecPolynomial{C}(a, MonomialVector{C}(x))

(::Type{VecPolynomial{C}}){C}(α) = VecPolynomial(Term{C}(α))
VecPolynomial{C}(x::PolyType{C}) = VecPolynomial{C}(x)
VecPolynomial(p::VecPolynomial) = p
VecPolynomial{C, T}(t::Term{C, T}) = VecPolynomial{C, T}([t.α], [t.x])
Base.convert{C, T}(::Type{VecPolynomial{C, T}}, x) = VecPolynomial(Term{C, T}(x))
Base.convert{C, T}(::Type{VecPolynomial{C, T}}, t::Term{C}) = VecPolynomial{C, T}([T(t.α)], [t.x])
Base.convert{C, T}(::Type{VecPolynomial{C, T}}, p::VecPolynomial{C, T}) = p
Base.convert{C, S, T}(::Type{VecPolynomial{C, T}}, p::VecPolynomial{C, S}) = VecPolynomial{C}(Vector{T}(p.a), p.x)

Base.convert{C, T}(::Type{TermContainer{C, T}}, p::VecPolynomial{C}) = VecPolynomial{C, T}(p)

function (::Type{VecPolynomial{C, T}}){C, T}(f::Function, x::MonomialVector{C})
    a = T[f(i) for i in 1:length(x)]
    VecPolynomial{C, T}(a, x)
end
(::Type{VecPolynomial{C, T}}){C, T}(f::Function, x::Vector) = VecPolynomial{C, T}(f, MonomialVector{C}(x))

Base.convert(::Type{Any}, p::VecPolynomial) = p

Base.convert{C}(::Type{PolyType{C}}, p::TermContainer{C}) = p
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

vars(p::VecPolynomial) = vars(p.x)

length(p::VecPolynomial) = length(p.a)
isempty(p::VecPolynomial) = length(p) > 0
start(::VecPolynomial) = 1
done(p::VecPolynomial, state) = length(p) < state
next(p::VecPolynomial, state) = (p[state], state+1)

extdeg(p::VecPolynomial) = extdeg(p.x)
mindeg(p::VecPolynomial) = mindeg(p.x)
maxdeg(p::VecPolynomial) = maxdeg(p.x)

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
    P(p.a[I], p.x[I])
end

function removedups{T}(adup::Vector{T}, Zdup::Vector{Vector{Int}})
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
    a, Z
end
function vecpolynomialclean{C, T}(vars::Vector{PolyVar{C}}, adup::Vector{T}, Zdup::Vector{Vector{Int}})
    a, Z = removedups(adup, Zdup)
    VecPolynomial{C, T}(a, MonomialVector{C}(vars, Z))
end

eltype{C, T}(::Type{VecPolynomial{C, T}}) = T
getindex(p::VecPolynomial, I::Int) = Term(p.a[I[1]], p.x[I[1]])
removemonomials{C}(p::VecPolynomial{C}, x::Vector) = removemonomials(p, MonomialVector{C}(x))

type MatPolynomial{C, T} <: TermType{C, T}
    Q::Vector{T}
    x::MonomialVector{C}
end

function trimap(i, j, n)
    div(n*(n+1), 2) - div((n-i+1)*(n-i+2), 2) + j-i+1
end

function trimat{T}(::Type{T}, f, n)
    Q = Vector{T}(trimap(n, n, n))
    for i in 1:n
        for j in i:n
            Q[trimap(i, j, n)] = f(i, j)
        end
    end
    Q
end

function getindex(p::MatPolynomial, I::NTuple{2,Int})
    i, j = I
    if i < j
        i, j = (j, i)
    end
    n = length(p.x)
    p.Q[trimap(i,j,n)]
end

function getmat{C, T}(p::MatPolynomial{C, T})
    n = length(p.x)
    A = Matrix{T}(n, n)
    for i in 1:n, j in i:n
        A[j,i] = A[i,j] = p.Q[trimap(i,j,n)]
    end
    A
end

function (::Type{MatPolynomial{C, T}}){C, T}(f::Function, x::MonomialVector{C})
    MatPolynomial{C, T}(trimat(T, f, length(x)), x)
end
(::Type{MatPolynomial{C, T}}){C, T}(f::Function, x::Vector) = MatPolynomial{C, T}(f, MonomialVector{C}(x))

function MatPolynomial{C, T}(Q::Matrix{T}, x::MonomialVector{C})
    MatPolynomial{C, T}((i,j) -> Q[i,j], x)
end
MatPolynomial{T<:VectorOfPolyType{false}}(Q::Matrix, x::Vector{T}) = MatPolynomial(Q, MonomialVector{false}(x))
MatPolynomial{T<:VectorOfPolyType{true}}(Q::Matrix, x::Vector{T}) = MatPolynomial(Q, MonomialVector{true}(x))

function Base.convert{C, T}(::Type{VecPolynomial{C, T}}, p::MatPolynomial{C, T})
    if isempty(p.Q)
        zero(VecPolynomial{C, T})
    else
        n = length(p.x)
        U = typeof(2*p.Q[1] + p.Q[1])
        if C
            N = trimap(n, n, n)
            Z = Vector{Vector{Int}}(N)
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
            v = p.x.vars
        else
            N = n^2
            x = Vector{Monomial}(N)
            a = Vector{U}(N)
            offset = 0
            for i in 1:n
                # for j in 1:n wouldn't be cache friendly for p.Q
                for j in i:n
                    k = trimap(i, j, n)
                    q = p.Q[k]
                    x[offset+k] = p.x[i] * p.x[j]
                    a[offset+k] = q
                    if i != j
                        offset += 1
                        x[offset+k] = p.x[j] * p.x[i]
                        a[offset+k] = q
                    end
                end
            end
            perm, X = sortmonovec(PolyVar{false}, x)
            a = a[perm]
            v = X.vars
            Z = X.Z
        end
        vecpolynomialclean(v, a, Z)
    end
end
VecPolynomial{C, T}(p::MatPolynomial{C, T}) = convert(VecPolynomial{C, T}, p)
TermContainer(p::MatPolynomial) = VecPolynomial(p)

type SOSDecomposition{C, T} <: TermType{C, T}
    ps::Vector{VecPolynomial{C, T}}
    function SOSDecomposition(ps::Vector{VecPolynomial{C, T}})
        new(ps)
    end
end
function (::Type{SOSDecomposition{C, T}}){C, T}(ps::Vector)
    SOSDecomposition(Vector{VecPolynomial{C, T}}(ps))
end
function (::Type{SOSDecomposition{C}}){C}(ps::Vector)
    T = reduce(promote_type, Int, map(eltype, ps))
    SOSDecomposition{C, T}(ps)
end
SOSDecomposition{T<:VectorOfPolyType{false}}(ps::Vector{T}) = SOSDecomposition{false}(ps)
SOSDecomposition{T<:VectorOfPolyType{true}}(ps::Vector{T}) = SOSDecomposition{true}(ps)

function Base.convert{C, T}(::Type{MatPolynomial{C, T}}, p::SOSDecomposition)
    X = mergemonovec(map(p -> p.x, p))
    m = length(p)
    n = length(X)
    Q = zeros(T, m, n)
    for i in 1:m
        q = p[i]
        k = 1
        for j in 1:n
            if k <= length(q) && X[j] == q.x[k]
                Q[i, j] = q.a[k]
                k += 1
            end
        end
    end
    MatPolynomial(Q' * Q, X)
end
MatPolynomial{C, T}(p::SOSDecomposition{C, T}) = MatPolynomial{C, T}(p)

function Base.convert{C, T}(::Type{SOSDecomposition{C, T}}, p::MatPolynomial)
    n = length(p.x)
    # TODO LDL^T factorization for SDP is missing in Julia
    # it would be nice to have though
    A = getmat(p)
    Q = chol(A)
    m = size(Q, 1)
    ps = [VecPolynomial{C, T}(Q[i,:], p.x) for i in 1:m]
    SOSDecomposition(ps)
end
# Without LDL^T, we need to do float(T)
SOSDecomposition{C, T}(p::MatPolynomial{C, T}) = SOSDecomposition{C, float(T)}(p)

length(p::SOSDecomposition) = length(p.ps)
isempty(p::SOSDecomposition) = isempty(p.ps)
start(p::SOSDecomposition) = start(p.ps)
done(p::SOSDecomposition, state) = done(p.ps, state)
next(p::SOSDecomposition, state) = next(p.ps, state)
getindex(p::SOSDecomposition, i::Int) = p.ps[i]

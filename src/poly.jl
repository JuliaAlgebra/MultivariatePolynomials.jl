export getmat, monomials, removemonomials
typealias AbstractTermContainer{T} Union{TermContainer{T}, NCTermContainer{T}}
typealias AbstractTermType{T} Union{TermType{T}, NCTermType{T}}
typealias AbstractTerm{T} Union{Term{T}, NCTerm{T}}

eltype{T}(p::AbstractTermType{T}) = T

Base.convert(::Type{Any}, t::AbstractTerm) = t
Base.copy{T<:AbstractTerm}(t::T) = T(copy(t.α), copy(t.x))

vars(t::AbstractTerm) = vars(t.x)

eltype{T}(::Type{AbstractTerm{T}}) = T
length(::AbstractTerm) = 1
isempty(::AbstractTerm) = false
start(::AbstractTerm) = false
done(::AbstractTerm, state) = state
next(t::AbstractTerm, state) = (t, true)
getindex(t::AbstractTerm, I::Int) = t

typealias AbstractVecPolynomial{T} Union{VecPolynomial{T}, NCVecPolynomial{T}}
Base.convert(::Type{Any}, p::AbstractVecPolynomial) = p

Base.convert(::Type{PolyType}, p::AbstractTermContainer) = p
function Base.convert{S}(::Type{S}, p::AbstractTermContainer)
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

vars(p::AbstractVecPolynomial) = vars(p.x)

length(p::AbstractVecPolynomial) = length(p.a)
isempty(p::AbstractVecPolynomial) = length(p) > 0
start(::AbstractVecPolynomial) = 1
done(p::AbstractVecPolynomial, state) = length(p) < state
next(p::AbstractVecPolynomial, state) = (p[state], state+1)

extdeg(p::AbstractVecPolynomial) = extdeg(p.x)
mindeg(p::AbstractVecPolynomial) = mindeg(p.x)
maxdeg(p::AbstractVecPolynomial) = maxdeg(p.x)

monomials(p::AbstractVecPolynomial) = p.x
function removemonomials(p::AbstractVecPolynomial, x::AbstractMonomialVector)
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

typealias AbstractMatPolynomial{T} Union{MatPolynomial{T}, NCMatPolynomial{T}}

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

function getindex(p::AbstractMatPolynomial, I::NTuple{2,Int})
    i, j = I
    if i < j
        i, j = (j, i)
    end
    n = length(p.x)
    p.Q[trimap(i,j,n)]
end

function getmat{T}(p::AbstractMatPolynomial{T})
    n = length(p.x)
    A = Matrix{T}(n, n)
    for i in 1:n, j in i:n
        A[j,i] = A[i,j] = p.Q[trimap(i,j,n)]
    end
    A
end

typealias AbstractSOSDecomposition{T} Union{SOSDecomposition{T}, NCSOSDecomposition{T}}

length(p::AbstractSOSDecomposition) = length(p.ps)
isempty(p::AbstractSOSDecomposition) = isempty(p.ps)
start(p::AbstractSOSDecomposition) = start(p.ps)
done(p::AbstractSOSDecomposition, state) = done(p.ps, state)
next(p::AbstractSOSDecomposition, state) = next(p.ps, state)

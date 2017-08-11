export MatPolynomial, SOSDecomposition
export getmat

struct SymMatrix{T} <: AbstractMatrix{T}
    Q::Vector{T}
    n
end

# i < j
function trimap(i, j, n)
    div(n*(n+1), 2) - div((n-i+1)*(n-i+2), 2) + j-i+1
end

function trimat(::Type{T}, f, n, σ) where {T}
    Q = Vector{T}(trimap(n, n, n))
    for i in 1:n
        for j in i:n
            Q[trimap(i, j, n)] = f(σ[i], σ[j])
        end
    end
    SymMatrix{T}(Q, n)
end

Base.size(Q::SymMatrix) = (Q.n, Q.n)

function Base.getindex(Q::SymMatrix, i, j)
    Q.Q[trimap(min(i, j), max(i, j), Q.n)]
end
function Base.getindex(Q::SymMatrix, k)
    i, j = divrem(k-1, Q.n)
    Q[i+1, j+1]
end
Base.getindex(Q::SymMatrix, I::Tuple) = Q[I...]

struct MatPolynomial{T, MT <: AbstractMonomial, MVT <: AbstractVector{MT}} <: AbstractPolynomialLike{T} # should be AbstractPolynomialLike{eltype(T)} but it doesn't work
    Q::SymMatrix{T}
    x::MVT
end
# When taking the promotion of a MatPolynomial of JuMP.Variable with a Polynomial JuMP.Variable, it should be a Polynomial of AffExpr
coefficienttype(::Type{<:MatPolynomial{T}}) where {T} = Base.promote_op(+, T, T)
polynomialtype(::Type{MatPolynomial{T, MT, MVT}}) where {T, MT, MVT} = polynomialtype(MT, coefficienttype(MatPolynomial{T, MT, MVT}))
polynomialtype(::Type{MatPolynomial{T, MT, MVT}}, ::Type{S}) where {S, T, MT, MVT} = polynomialtype(MT, S)

Base.zero(::Type{MatPolynomial{T, MT, MVT}}) where {T, MT, MVT} = MatPolynomial{T, MT, monovectype(MT)}(SymMatrix{T}(T[], 0), emptymonovec(MT))

Base.getindex(p::MatPolynomial, I...) = getindex(p.Q, I...)

getmat(p::MatPolynomial{T}) where {T} = p.Q

function MatPolynomial{T}(f::Function, x::AbstractVector{MT}, σ) where {T, MT<:AbstractMonomial}
    MatPolynomial{T, MT, monovectype(x)}(trimat(T, f, length(x), σ), x)
end
MatPolynomial{T}(f::Function, x::AbstractVector, σ) where T = MatPolynomial{T}(f, monomial.(x), σ)
function MatPolynomial{T}(f::Function, x::AbstractVector) where T
    σ, X = sortmonovec(x)
    MatPolynomial{T}(f, X, σ)
end
MatPolynomial(f::Function, x) = MatPolynomial{Base.promote_op(f, Int, Int)}(f, x)

function MatPolynomial(Q::AbstractMatrix{T}, x, σ) where {T}
    MatPolynomial{T}((i,j) -> Q[σ[i], σ[j]], x)
end
function MatPolynomial(Q::AbstractMatrix{T}, x) where {T}
    σ, X = sortmonovec(x)
    MatPolynomial(Q, X, σ)
end

#function Base.convert{T, PT <: AbstractPolynomial{T}}(::Type{PT}, p::MatPolynomial)
#    # coefficienttype(p) may be different than T and polynomial(p) may be different than PT (different module)
#    convert(PT, polynomial(p))
#end
function polynomial(p::MatPolynomial)
    polynomial(getmat(p), p.x)
end
function polynomial(p::MatPolynomial, ::Type{S}) where {S}
    polynomial(getmat(p), p.x, S)
end

struct SOSDecomposition{T, PT <: AbstractPolynomialLike{T}} <: AbstractPolynomialLike{T} # If promote_op((x, y) -> x * y + x * y, T, T) != T then it might not be true
    ps::Vector{PT}
    function SOSDecomposition{T, PT}(ps::Vector{PT}) where {T, PT}
        new(ps)
    end
end
SOSDecomposition(ps::Vector{PT}) where {T, PT <: AbstractPolynomialLike{T}} = SOSDecomposition{T, PT}(ps)
polynomialtype(::Type{SOSDecomposition{T, PT}}) where {T, PT} = polynomialtype(PT)
#function SOSDecomposition(ps::Vector)
#    T = reduce(promote_type, Int, map(eltype, ps))
#    SOSDecomposition{T}(ps)
#end

function MatPolynomial(p::SOSDecomposition{T}) where {T}
    X = mergemonovec(map(monomials, p))
    m = length(p)
    n = length(X)
    Q = zeros(T, m, n)
    for i in 1:m
        j = 1
        for t in terms(p[i])
            while X[j] != monomial(t)
                j += 1
            end
            Q[i, j] = coefficient(t)
            j += 1
        end
    end
    MatPolynomial(Q' * Q, X)
end

function SOSDecomposition(p::MatPolynomial)
    n = length(p.x)
    # TODO LDL^T factorization for SDP is missing in Julia
    # it would be nice to have though
    A = getmat(p)
    Q = chol(A)
    m = size(Q, 1)
    ps = [polynomial(Q[i,:], p.x) for i in 1:m]
    SOSDecomposition(ps)
end
# Without LDL^T, we need to do float(T)
SOSDecomposition(p::MatPolynomial{C, T}) where {C, T} = SOSDecomposition{C, float(T)}(p)

Base.length(p::SOSDecomposition) = length(p.ps)
Base.isempty(p::SOSDecomposition) = isempty(p.ps)
Base.start(p::SOSDecomposition) = start(p.ps)
Base.done(p::SOSDecomposition, state) = done(p.ps, state)
Base.next(p::SOSDecomposition, state) = next(p.ps, state)
Base.getindex(p::SOSDecomposition, i::Int) = p.ps[i]

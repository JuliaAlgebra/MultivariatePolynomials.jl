export Moment, Measure, MatMeasure, AtomicMeasure, zeta, ζ, extractatoms
using RowEchelon

type Moment{C, T}
    α::T
    x::Monomial{C}
end

# If a monomial is not in x, it does not mean that the moment is zero, it means that it is unknown/undefined
type Measure{C, T}
    a::Vector{T}
    x::MonomialVector{C}

    function Measure(a::Vector{T}, x::MonomialVector{C})
        if length(a) != length(x)
            error("There should be as many coefficient than monomials")
        end
        new(a, x)
    end
end

Measure{C, T}(a::Vector{T}, x::MonomialVector{C}) = Measure{C, T}(a, x)
function (::Type{Measure{C}}){C}(a::Vector, x::Vector)
    if length(a) != length(x)
        error("There should be as many coefficient than monomials")
    end
    σ, X = sortmonovec(PolyVar{C}, x)
    Measure(a[σ], X)
end
Measure{T<:VectorOfPolyType{true}}(a::Vector, X::Vector{T}) = Measure{true}(a, X)
Measure{T<:VectorOfPolyType{false}}(a::Vector, X::Vector{T}) = Measure{false}(a, X)

function (*)(α, μ::Measure)
    Measure(α * μ.a, μ.x)
end
function (+)(μ::Measure, ν::Measure)
    if μ.x != ν.x
        # Should I just drop moments that are not common ?
        throw(ArgumentError("Cannot sum measures with different moments"))
    end
    Measure(μ.a + ν.a, μ.x)
end

function ζ{C, T}(v::Vector{T}, x::MonomialVector{C}, varorder::Vector{PolyVar{C}})
    Measure(T[m(v, varorder) for m in x], x)
end

type MatMeasure{C, T}
    Q::Vector{T}
    x::MonomialVector{C}
end
vars{M<:Union{Moment, Measure, MatMeasure}}(μ::M) = vars(μ.x)
nvars{M<:Union{Moment, Measure, MatMeasure}}(μ::M) = nvars(μ.x)
function (::Type{MatMeasure{C, T}}){C, T}(f::Function, x::MonomialVector{C}, σ=1:length(x))
    MatMeasure{C, T}(trimat(T, f, length(x), σ), x)
end
function (::Type{MatMeasure{C, T}}){C, T}(f::Function, x::Vector)
    σ, X = sortmonovec(x)
    MatMeasure{C, T}(f, X, σ)
end
(::Type{MatMeasure{C}}){C}(f::Function, x) = MatMeasure{C, Base.promote_op(f, Int, Int)}(f, x)
MatMeasure{T<:VectorOfPolyType{false}}(f::Function, x::Vector{T}) = MatMeasure{false}(f, x)
MatMeasure{T<:VectorOfPolyType{true}}(f::Function, x::Vector{T}) = MatMeasure{true}(f, x)
MatMeasure{C}(f::Function, x::MonomialVector{C}) = MatMeasure{C}(f, x)

MatMeasure{C, T}(μ::Measure{C, T}, x) = MatMeasure{C, T}(μ, x)
(::Type{MatMeasure{C, T}}){C, T}(μ::Measure{C, T}, x) = MatMeasure{C, T}(μ, MonomialVector(x))
function (::Type{MatMeasure{C, T}}){C, T}(μ::Measure{C, T}, x::MonomialVector{C})
    function getmom(i, j)
        k = multisearch(μ.x, MonomialVector([x[i]*x[j]]))[1]
        if k == 0
            throw(ArgumentError("μ does not have the moment $(x[i]*x[j])"))
        end
        μ.a[k]
    end
    MatMeasure{C, T}(getmom, x)
end

function MatMeasure{C, T}(Q::Matrix{T}, x::MonomialVector{C})
    MatMeasure{C, T}((i,j) -> Q[i, j], x)
end
function matmeasperm{C, T}(Q::Matrix{T}, x::MonomialVector{C}, σ)
    MatMeasure{C, T}((i,j) -> Q[σ[i], σ[j]], x)
end
function MatMeasure{T}(Q::Matrix{T}, x::Vector)
    σ, X = sortmonovec(x)
    matmeasperm(Q, X, σ)
end

function getmat{C, T}(μ::MatMeasure{C, T})
    _getmat(μ.Q, length(μ.x))
end

type AtomicMeasure{C, T}
    vars::Vector{PolyVar{C}}
    λ::Vector{T}
    vals::Vector{Vector{T}}
end
function AtomicMeasure{C, S, T}(vars::Vector{PolyVar{C}}, λ::Vector{S}, vals::Vector{Vector{T}})
    AtomicMeasure{C, promote_type(S, T)}(vars, λ, vals)
end

Measure{C, T}(μ::AtomicMeasure{C, T}, x) = Measure{C, T}(μ, x)
(::Type{Measure{C, T}}){C, T}(μ::AtomicMeasure{C, T}, x::Vector) = Measure{C, T}(μ, MonomialVector(x))
function (::Type{Measure{C, T}}){C, T}(μ::AtomicMeasure{C, T}, x::MonomialVector)
    sum(μ.λ[i] * ζ(μ.vals[i], x, μ.vars) for i in 1:length(μ.vals))
end

# Solve the system
# y = [U 0] * y
# where y = x[end:-1:1]
# which is
# y = U * β
function solve_system(U, x)
    m, r = size(U)
    @assert m == length(x)
    n = nvars(x)
    v = vars(x)
    pivots = [findfirst(j -> U[j, i] != 0, 1:m) for i in 1:r]
    if any(pivots .== 0)
        keep = pivots .> 0
        pivots = pivots[keep]
        r = length(pivots)
        U = U[:, keep]
    end
    β = x[m+1-pivots]
    function multisearch_check(y)
        idxs = multisearch(x, y)
        if any(idxs .== 0)
            error("Missing monomials $(y[idxs .== 0]) in $(x)")
        end
        idxs
    end
    Ns = [U[m+1-reverse(multisearch_check(v[i] * β)), :] for i in 1:n]
    λ = rand(n)
    λ /= sum(λ)
    N = sum(λ .* Ns)
    Z = schurfact(N)[:Z]
    vals = [Vector{Float64}(n) for j in 1:r]
    for j in 1:r
        qj = Z[:, j]
        for i in 1:n
            vals[j][i] = dot(qj, Ns[i] * qj)
        end
    end
    r, vals
end

function extractatoms(μ::MatMeasure, tol::Real, shift::Real, ɛ::Real=-1)
    # We reverse the ordering so that the first columns corresponds to low order monomials
    # so that we have more chance that low order monomials are in β and then more chance
    # v[i] * β to be in μ.x
    M = getmat(μ)[end:-1:1, end:-1:1]
    m = size(M, 1)
    #r = rank(M, tol)
    U = chol(M + shift * eye(m))
    V = U[map(i -> abs(U[i, i]) > tol, 1:m), :]
    r = size(V, 1)
#   F = svdfact(M)
#   S = F.S
#   r = sum(F.S .> tol)
#   V = F.U[:, 1:r] .* repmat(sqrt.(S[1:r])', size(F.U, 1), 1)
    rref!(V, ɛ == -1 ? sqrt(eps(norm(V, Inf))) : ɛ)
    r, vals = solve_system(V', μ.x)
    # Determine weights
    Ms = similar(M, r, r)
    v = vars(μ)
    for i in 1:r
        vi = ζ(vals[i], μ.x, v)
        Ms[:, i] = vi.a[end:-1:end-r+1] * vi.a[end]
    end
    λ = Ms \ M[1:r, 1]
    AtomicMeasure(v, λ, vals)
end

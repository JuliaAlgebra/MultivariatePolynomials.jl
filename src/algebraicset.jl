export AbstractSemialgebraicSet, AbstractBasicSemialgebraicSet, AbstractAlgebraicSet
export FullSpace, AlgebraicSet, BasicSemialgebraicSet, addequality!, addinequality!
# Semialgebraic set described by polynomials with coefficients in T
abstract type AbstractSemialgebraicSet end

abstract type AbstractBasicSemialgebraicSet <: AbstractSemialgebraicSet end
abstract type AbstractAlgebraicSet <: AbstractBasicSemialgebraicSet end

addinequality!(S::AbstractAlgebraicSet, p) = throw(ArgumentError("Cannot add inequality to an algebraic set"))

immutable FullSpace <: AbstractAlgebraicSet
end

type AlgebraicSet{T, PT<:APL{T}} <: AbstractAlgebraicSet
    p::Vector{PT}
end
function AlgebraicSet{T, PT}() where {T, PT<:APL{T}}
    AlgebraicSet(PT[])
end

addequality!(V::AlgebraicSet, p) = push!(V.p, p)
Base.intersect(S::AlgebraicSet, T::AlgebraicSet) = AlgebraicSet([S.p; T.p])

type BasicSemialgebraicSet{T, PT<:APL{T}} <: AbstractBasicSemialgebraicSet
    V::AlgebraicSet{T, PT}
    p::Vector{PT}
end
function BasicSemialgebraicSet{T, PT}() where {T, PT<:APL{T}}
    BasicSemialgebraicSet(AlgebraicSet{T, PT}(), PT[])
end

addequality!(S::BasicSemialgebraicSet, p) = addequality!(S.V, p)
addinequality!(S::BasicSemialgebraicSet, p) = push!(S.p, p)

Base.intersect(S::BasicSemialgebraicSet, T::BasicSemialgebraicSet) = BasicSemialgebraicSet(S.V ∩ T.V, [S.p; T.p])
Base.intersect(S::BasicSemialgebraicSet, T::AlgebraicSet) = BasicSemialgebraicSet(S.V ∩ T, copy(S.p))
Base.intersect(T::AlgebraicSet, S::BasicSemialgebraicSet) = intersect(S, T)

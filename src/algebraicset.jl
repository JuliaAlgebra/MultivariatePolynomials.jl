export AbstractSemialgebraicSet, AbstractBasicSemialgebraicSet, AbstractAlgebraicSet
export FullSpace, AlgebraicSet, BasicSemialgebraicSet, addequality!, addinequality!
# Semialgebraic set described by polynomials with coefficients in T
@compat abstract type AbstractSemialgebraicSet end

@compat abstract type AbstractBasicSemialgebraicSet <: AbstractSemialgebraicSet end
@compat abstract type AbstractAlgebraicSet <: AbstractBasicSemialgebraicSet end

addinequality!(S::AbstractAlgebraicSet, p) = throw(ArgumentError("Cannot add inequality to an algebraic set"))

immutable FullSpace <: AbstractAlgebraicSet
end

type AlgebraicSet <: AbstractAlgebraicSet
    p::Vector
end
function (::Type{AlgebraicSet})()
    AlgebraicSet(Any[])
end

addequality!(V::AlgebraicSet, p) = push!(V.p, p)
Base.intersect(S::AlgebraicSet, T::AlgebraicSet) = AlgebraicSet([S.p; T.p])

type BasicSemialgebraicSet <: AbstractBasicSemialgebraicSet
    V::AlgebraicSet
    p::Vector
end
function (::Type{BasicSemialgebraicSet})()
    BasicSemialgebraicSet(AlgebraicSet(), Any[])
end

addequality!(S::BasicSemialgebraicSet, p) = addequality!(S.V, p)
addinequality!(S::BasicSemialgebraicSet, p) = push!(S.p, p)

Base.intersect(S::BasicSemialgebraicSet, T::BasicSemialgebraicSet) = BasicSemialgebraicSet(S.V ∩ T.V, [S.p; T.p])
Base.intersect(S::BasicSemialgebraicSet, T::AlgebraicSet) = BasicSemialgebraicSet(S.V ∩ T, copy(S.p))
Base.intersect(T::AlgebraicSet, S::BasicSemialgebraicSet) = intersect(S, T)

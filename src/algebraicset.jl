export AlgebraicSet, BasicSemialgebraicSet, addequality!, addinequality!
# Semialgebraic set described by polynomials with coefficients in T
abstract AbstractSemialgebraicSet

abstract AbstractBasicSemialgebraicSet <: AbstractSemialgebraicSet
abstract AbstractAlgebraicSet <: AbstractBasicSemialgebraicSet

addinequality!(S::AbstractAlgebraicSet, p) = throw(ArgumentError("Cannot add inequality to an algebraic set"))

type AlgebraicSet <: AbstractAlgebraicSet
    p::Vector
end
function (::Type{AlgebraicSet})()
    AlgebraicSet(Any[])
end

addequality!(V::AlgebraicSet, p) = push!(V.p, p)

type BasicSemialgebraicSet <: AbstractBasicSemialgebraicSet
    V::AlgebraicSet
    p::Vector
end
function (::Type{BasicSemialgebraicSet})()
    BasicSemialgebraicSet(AlgebraicSet(), Any[])
end

addequality!(S::BasicSemialgebraicSet, p) = addequality!(S.V, p)
addinequality!(S::BasicSemialgebraicSet, p) = push!(S.p, p)

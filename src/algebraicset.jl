export AlgebraicSet, BasicSemialgebraicSet, addequality!, addinequality!
# Semialgebraic set described by polynomials with coefficients in T
abstract AbstractSemialgebraicSet{T}

addequality!{T}(S::AbstractSemialgebraicSet{T}, p) = addequality!(S, VecPolynomial{T}(p))
addinequality!{T}(S::AbstractSemialgebraicSet{T}, p) = addinequality!(S, VecPolynomial{T}(p))

abstract AbstractBasicSemialgebraicSet{T} <: AbstractSemialgebraicSet{T}
abstract AbstractAlgebraicSet{T} <: AbstractBasicSemialgebraicSet{T}

addinequality!{T}(S::AbstractAlgebraicSet{T}, p) = throw(ArgumentError("Cannot add inequality to an algebraic set"))

type AlgebraicSet{T} <: AbstractAlgebraicSet{T}
    p::Vector{VecPolynomial{T}}
end
function (::Type{AlgebraicSet{T}}){T}()
    AlgebraicSet{T}(VecPolynomial{T}[])
end


function Base.convert{T,U}(::Type{AlgebraicSet{T}}, V::AlgebraicSet{U})
    AlgebraicSet{T}(Vector{VecPolynomial{T}}(V.p))
end

addequality!{T}(V::AlgebraicSet{T}, p::VecPolynomial{T}) = push!(V.p, p)

type BasicSemialgebraicSet{T} <: AbstractBasicSemialgebraicSet{T}
    V::AlgebraicSet{T}
    p::Vector{VecPolynomial{T}}
end
function (::Type{BasicSemialgebraicSet{T}}){T}()
    BasicSemialgebraicSet{T}(AlgebraicSet{T}(), VecPolynomial{T}[])
end


function Base.convert{T,U}(::Type{BasicSemialgebraicSet{T}}, S::BasicSemialgebraicSet{U})
    BasicSemialgebraicSet{T}(convert(AlgebraicSet{T}, S.V), Vector{VecPolynomial{T}}(S.p))
end

addequality!{T}(S::BasicSemialgebraicSet{T}, p::VecPolynomial{T}) = addequality!(S.V, p)
addinequality!{T}(S::BasicSemialgebraicSet{T}, p::VecPolynomial{T}) = push!(S.p, p)

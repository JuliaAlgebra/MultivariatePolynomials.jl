export RationalPoly
import Base.+, Base.-, Base.*, Base./

immutable RationalPoly{S,T} <: PolyType
    num::TermContainer{S}
    den::TermContainer{T}
end

Base.convert{S,T}(::Type{RationalPoly{S,T}}, q::RationalPoly{S,T}) = q
Base.convert{S,T,U,V}(::Type{RationalPoly{S,T}}, q::RationalPoly{U,V}) = TermContainer{S}(q.num) / TermContainer{T}(q.den)
function Base.convert{S,T}(::Type{RationalPoly{S,T}}, p::TermContainer{S})
    p / one(TermContainer{T})
end
function Base.convert{S,T}(::Type{RationalPoly{S,T}}, p)
    Base.convert(RationalPoly{S,T}, TermContainer{S}(p))
end

function (/){S,T}(num::TermContainer{S}, den::TermContainer{T})
    RationalPoly{S,T}(num, den)
end
(/)(num::PolyType, den) = num / TermContainer(den)
function (/)(num, den::PolyType)
    TermContainer(num) / den
end
(/)(num::PolyType, den::PolyType) = TermContainer(num) / TermContainer(den)

function (+)(r::RationalPoly, s::RationalPoly)
    (r.num*s.den + r.den*s.num) / (r.den * s.den)
end
function (+)(p::VecPolynomial, r::RationalPoly)
    (p*r.den + r.num) / r.den
end
(+)(r::RationalPoly, p::VecPolynomial) = p + r
function (-)(r::RationalPoly, s::RationalPoly)
    (r.num*s.den - r.den*s.num) / (r.den * s.den)
end
(-)(p::PolyType, s::RationalPoly) = (p * s.den - s.num) / s.den
(-)(s::RationalPoly, p::PolyType) = (s.num - p * s.den) / s.den

(*)(r::RationalPoly, s::RationalPoly) = (r.num*s.num) / (r.den*s.den)
function (*)(p::TermContainer, r::RationalPoly)
    if p == r.den
        r.num
    else
        (p * r.num) / r.den
    end
end
(*)(r::RationalPoly, p::Term) = p * r
(*)(r::RationalPoly, p::VecPolynomial) = p * r
(*)(p::PolyType, r::RationalPoly) = TermContainer(p) * r
(*)(r::RationalPoly, p::PolyType) = r * TermContainer(p)
(*)(α, r::RationalPoly) = TermContainer(α) * TermContainer(p)
(*)(r::RationalPoly, α) = r * TermContainer(α)

zero(r::RationalPoly) = zero(r.num)
zero{T<:RationalPoly}(::Type{T}) = zero(VecPolynomial)
zero{T}(::Type{RationalPoly{T}}) = zero(VecPolynomial{T})

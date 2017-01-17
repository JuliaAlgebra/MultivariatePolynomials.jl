export RationalPoly
import Base.+, Base.-, Base.*, Base./

immutable RationalPoly{C, S, T} <: PolyType{C}
    num::TermContainer{C, S}
    den::TermContainer{C, T}
end
iscomm{C, S, T}(r::Type{RationalPoly{C, S, T}}) = C

Base.convert{C, S, T}(::Type{RationalPoly{C, S, T}}, q::RationalPoly{C, S, T}) = q
Base.convert{C, S, T, U, V}(::Type{RationalPoly{C, S, T}}, q::RationalPoly{C, U, V}) = TermContainer{C, S}(q.num) / TermContainer{C, T}(q.den)
function Base.convert{C, S, T}(::Type{RationalPoly{C, S, T}}, p::TermContainer{C, S})
    p / one(TermContainer{C, T})
end
function Base.convert{C, S, T}(::Type{RationalPoly{C, S, T}}, p)
    Base.convert(RationalPoly{C, S, T}, TermContainer{C, S}(p))
end

function (/){C, S, T}(num::TermContainer{C, S}, den::TermContainer{C, T})
    RationalPoly{C, S, T}(num, den)
end
(/){C}(num::PolyType{C}, den) = num / TermContainer{C}(den)
function (/){C}(num, den::PolyType{C})
    TermContainer{C}(num) / den
end
(/){C}(num::PolyType{C}, den::PolyType{C}) = TermContainer{C}(num) / TermContainer{C}(den)

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
(*){C}(α, r::RationalPoly{C}) = TermContainer{C}(α) * TermContainer{C}(p)
(*){C}(r::RationalPoly{C}, α) = r * TermContainer{C}(α)

zero(r::RationalPoly) = zero(r.num)
zero{T<:RationalPoly}(::Type{T}) = zero(VecPolynomial)
zero{C, T}(::Type{RationalPoly{C, T}}) = zero(VecPolynomial{C, T})

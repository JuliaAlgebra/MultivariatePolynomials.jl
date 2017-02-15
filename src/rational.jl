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
function Base.convert{C, S, T}(::Type{RationalPoly{C, S, T}}, p::TermContainer)
    convert(RationalPoly{C, S, T}, TermContainer{C, S}(p))
end
function Base.convert{C, S, T}(::Type{RationalPoly{C, S, T}}, p)
    Base.convert(RationalPoly{C, S, T}, TermContainer{C, S}(p))
end

(/)(r::RationalPoly, p::TermContainer) = r.num / (r.den * p)
function (/){C, S, T}(num::TermContainer{C, S}, den::TermContainer{C, T})
    RationalPoly{C, S, T}(num, den)
end
function (/){C}(num, den::PolyType{C})
    TermContainer{C}(num) / den
end
(/){C}(num::PolyType{C}, den::PolyType{C}) = TermContainer{C}(num) / TermContainer{C}(den)

# Polynomial divided by coefficient is a polynomial not a rational polynomial
(/){C}(num::PolyType{C}, den) = num * (1 / den)

function (+)(r::RationalPoly, s::RationalPoly)
    (r.num*s.den + r.den*s.num) / (r.den * s.den)
end
function (+)(p::TermContainer, r::RationalPoly)
    (p*r.den + r.num) / r.den
end
(+)(r::RationalPoly, p::Polynomial) = p + r
(+)(r::RationalPoly, t::Term) = t + r
function (-)(r::RationalPoly, s::RationalPoly)
    (r.num*s.den - r.den*s.num) / (r.den * s.den)
end
(-)(p::PolyType, s::RationalPoly) = (p * s.den - s.num) / s.den
(-)(s::RationalPoly, p::PolyType) = (s.num - p * s.den) / s.den

(*)(r::RationalPoly, s::RationalPoly) = (r.num*s.num) / (r.den*s.den)
(*)(p::TermContainer, r::RationalPoly) = p == r.den ? r.num : (p * r.num) / r.den
(*)(r::RationalPoly, p::Polynomial) = p == r.den ? r.num : (r.num * p) / r.den
(*)(r::RationalPoly, t::Term)          = t == r.den ? r.num : (r.num * t) / r.den
(*)(p::PolyType, r::RationalPoly) = TermContainer(p) * r
(*)(r::RationalPoly, p::Monomial) = r * TermContainer(p)
(*)(r::RationalPoly, p::PolyVar)  = r * TermContainer(p)
(*){C}(α, r::RationalPoly{C}) = TermContainer{C}(α) * r
(*){C}(r::RationalPoly{C}, α) = r * TermContainer{C}(α)

zero(r::RationalPoly) = zero(r.num)
zero{C, S, T}(::Type{RationalPoly{C, S, T}}) = zero(Polynomial{C, S})

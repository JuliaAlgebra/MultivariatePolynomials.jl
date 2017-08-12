export RationalPoly
import Base.+, Base.-, Base.*, Base./

# TODO We should take gcd between numerator and denominator
struct RationalPoly{NT <: APL, DT <: APL}
    num::NT
    den::DT
end

Base.convert(::Type{RationalPoly{NT, DT}}, q::RationalPoly{NT, DT}) where {NT, DT} = q
Base.convert(::Type{RationalPoly{NTout, DTout}}, q::RationalPoly{NTin, DTin}) where {NTout, DTout, NTin, DTin} = convert(NTout, q.num) / convert(DTout, q.den)
#function Base.convert(::Type{RationalPoly{NT, DT}}, p::NT) where {NT, DT}
#    p / one(DT)
#end
function Base.convert(::Type{RationalPoly{NT, DT}}, p::APL) where {NT, DT}
    convert(NT, p) / one(DT)
end
function Base.convert(::Type{RationalPoly{NT, DT}}, α) where {NT, DT}
    convert(NT, α) / one(DT)
    #convert(RationalPoly{NT, DT}, convert(NT, α))
end

(/)(r::RationalPoly, p) = r.num / (r.den * p)
function (/)(num::NT, den::DT) where {NT <: APL, DT <: APL}
    RationalPoly{NT, DT}(num, den)
end
function (/)(num, den::APL)
    constantterm(num, den) / den
end
# Polynomial divided by coefficient is a polynomial not a rational polynomial
# (1/den) * num would not be correct in case of noncommutative coefficients
(/)(num::APL, den) = num * (1 / den)

function (+)(r::RationalPoly, s::RationalPoly)
    (r.num*s.den + r.den*s.num) / (r.den * s.den)
end
function _plus(r::RationalPoly, p)
    (p*r.den + r.num) / r.den
end
(+)(p::APL, r::RationalPoly) = _plus(r, p)
(+)(r::RationalPoly, p::APL) = _plus(r, p)
(+)(r::RationalPoly, α) = _plus(r, α)
(+)(α, r::RationalPoly) = _plus(r, α)
function (-)(r::RationalPoly, s::RationalPoly)
    (r.num*s.den - r.den*s.num) / (r.den * s.den)
end
_minus(p, s::RationalPoly) = (p * s.den - s.num) / s.den
_minus(s::RationalPoly, p) = (s.num - p * s.den) / s.den
(-)(p::APL, r::RationalPoly) = _minus(r, p)
(-)(r::RationalPoly, p::APL) = _minus(r, p)
(-)(r::RationalPoly, α) = _minus(r, α)
(-)(α, r::RationalPoly) = _minus(r, α)
(-)(r::RationalPoly) = (-r.num) / r.den

(*)(r::RationalPoly, s::RationalPoly) = (r.num*s.num) / (r.den*s.den)
(*)(p::APL, r::RationalPoly)          = (p * r.num) / r.den
(*)(r::RationalPoly, p::APL)          = (r.num * p) / r.den
(*)(α, r::RationalPoly)               = (α * r.num) / r.den
(*)(r::RationalPoly, α)               = (r.num * α) / r.den

Base.zero(::RationalPoly{NT}) where {NT} = zero(NT)
Base.zero(::Type{RationalPoly{NT, DT}}) where {NT, DT} = zero(NT)

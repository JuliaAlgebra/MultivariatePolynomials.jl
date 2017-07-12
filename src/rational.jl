export RationalPoly
import Base.+, Base.-, Base.*, Base./

immutable RationalPoly{NT <: APL, DT <: APL}
    num::NT
    den::DT
end

Base.convert{NT, DT}(::Type{RationalPoly{NT, DT}}, q::RationalPoly{NT, DT}) = q
Base.convert{NTout, DTout, NTin, DTin}(::Type{RationalPoly{NTout, DTout}}, q::RationalPoly{NTin, DTin}) = convert(NTout, q.num) / convert(DTout, q.den)
function Base.convert{NT, DT}(::Type{RationalPoly{NT, DT}}, p::NT)
    p / one(DT)
end
function Base.convert{NT, DT}(::Type{RationalPoly{NT, DT}}, p)
    convert(RationalPoly{NT, DT}, convert(NT, p))
end

(/)(r::RationalPoly, p) = r.num / (r.den * p)
function (/){NT <: APL, DT <: APL}(num::NT, den::DT)
    RationalPoly{NT, DT}(num, den)
end
function (/)(num, den::APL)
    term(num, den) / den
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
_minus(p::APL, s::RationalPoly) = (p * s.den - s.num) / s.den
_minus(s::RationalPoly, p::APL) = (s.num - p * s.den) / s.den
(-)(p::APL, r::RationalPoly) = _minus(r, p)
(-)(r::RationalPoly, p::APL) = _minus(r, p)
(-)(r::RationalPoly, α) = _minus(r, α)
(-)(α, r::RationalPoly) = _minus(r, α)

(*)(r::RationalPoly, s::RationalPoly) = (r.num*s.num) / (r.den*s.den)
# Not type stable, currently it is a hack for SumOfSquares/test/SOSdemo2.jl:line 22
# We should take gcd between numerator and denominator instead and in sosdemo2, we should cast to polynomial manually
(*)(p::APL, r::RationalPoly)          = p == r.den ? r.num : (p * r.num) / r.den
(*)(r::RationalPoly, p::APL)          = p == r.den ? r.num : (r.num * p) / r.den
(*)(α, r::RationalPoly)               = α == r.den ? r.num : (α * r.num) / r.den
(*)(r::RationalPoly, α)               = α == r.den ? r.num : (r.num * α) / r.den

Base.zero{NT}(::RationalPoly{NT}) = zero(NT)
Base.zero{NT, DT}(::Type{RationalPoly{NT, DT}}) = zero(NT)

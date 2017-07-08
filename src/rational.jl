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
    term(num) / den
end
# Polynomial divided by coefficient is a polynomial not a rational polynomial
# (1/den) * num would not be correct in case of noncommutative coefficients
(/)(num::APL, den) = num * (1 / den)

function (+)(r::RationalPoly, s::RationalPoly)
    (r.num*s.den + r.den*s.num) / (r.den * s.den)
end
function (+)(p, r::RationalPoly)
    (p*r.den + r.num) / r.den
end
(+)(r::RationalPoly, p) = p + r
function (-)(r::RationalPoly, s::RationalPoly)
    (r.num*s.den - r.den*s.num) / (r.den * s.den)
end
(-)(p, s::RationalPoly) = (p * s.den - s.num) / s.den
(-)(s::RationalPoly, p) = (s.num - p * s.den) / s.den

(*)(r::RationalPoly, s::RationalPoly) = (r.num*s.num) / (r.den*s.den)
# Not type stable, currently it is a hack for SumOfSquares/test/SOSdemo2.jl:line 22
# We should take gcd between numerator and denominator instead and in sosdemo2, we should cast to polynomial manually
(*)(p, r::RationalPoly)       = p == r.den ? r.num : (p * r.num) / r.den
(*)(r::RationalPoly, p)       = p == r.den ? r.num : (r.num * p) / r.den

zero{NT}(::RationalPoly{NT}) = zero(NT)
zero{NT, DT}(::Type{RationalPoly{NT, DT}}) = zero(NT)

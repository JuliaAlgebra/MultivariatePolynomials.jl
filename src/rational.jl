export RationalPoly

# TODO We should take gcd between numerator and denominator
struct RationalPoly{NT <: APL, DT <: APL}
    num::NT
    den::DT
end

# This constructor is called from LinearAlgebra in the method Matrix{T}(s::UniformScaling{Bool}, dims)
RationalPoly{NT,DT}(x::Bool) where {NT,DT} = ifelse(x, one(RationalPoly{NT,DT}), zero(RationalPoly{NT,DT}))

Base.numerator(r::RationalPoly) = r.num
Base.denominator(r::RationalPoly) = r.den

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

# This is called by the default implementation of `Base.oneunit`, still in Julia v1.8 at least
RationalPoly{NT,DT}(r::RationalPoly{NT,DT}) where {NT,DT} = r

Base.inv(r::RationalPoly) = r.den / r.num
Base.inv(p::APL{T}) where T = one(T) / p
Base.:/(r::RationalPoly, p) = r.num / (r.den * p)
Base.:/(r::RationalPoly, p::APL) = r.num / (r.den * p)
Base.:/(r::RationalPoly, s::RationalPoly) = (r.num * s.den) / (s.num * r.den)
function Base.:/(num::NT, den::DT) where {NT <: APL, DT <: APL}
    RationalPoly{NT, DT}(num, den)
end
function Base.:/(num, den::APL)
    constant_term(num, den) / den
end
# Polynomial divided by coefficient is a polynomial not a rational polynomial
# (1/den) * num would not be correct in case of noncommutative coefficients
Base.:/(num::APL, den) = map_coefficientsnz(α -> α/den, num)

function Base.:+(r::RationalPoly, s::RationalPoly)
    (r.num*s.den + r.den*s.num) / (r.den * s.den)
end
function _plus(r::RationalPoly, p)
    (p*r.den + r.num) / r.den
end
Base.:+(p::APL, r::RationalPoly) = _plus(r, p)
Base.:+(r::RationalPoly, p::APL) = _plus(r, p)
Base.:+(r::RationalPoly, α) = _plus(r, α)
Base.:+(α, r::RationalPoly) = _plus(r, α)
function Base.:-(r::RationalPoly, s::RationalPoly)
    (r.num*s.den - r.den*s.num) / (r.den * s.den)
end
_minus(p, s::RationalPoly) = (p * s.den - s.num) / s.den
_minus(s::RationalPoly, p) = (s.num - p * s.den) / s.den
Base.:-(p::APL, r::RationalPoly) = _minus(p, r)
Base.:-(r::RationalPoly, p::APL) = _minus(r, p)
Base.:-(r::RationalPoly, α) = _minus(r, α)
Base.:-(α, r::RationalPoly) = _minus(α, r)
Base.:-(r::RationalPoly) = (-r.num) / r.den

Base.:*(r::RationalPoly, s::RationalPoly) = (r.num*s.num) / (r.den*s.den)
Base.:*(p::APL, r::RationalPoly)          = (p * r.num) / r.den
Base.:*(r::RationalPoly, p::APL)          = (r.num * p) / r.den
Base.:*(α, r::RationalPoly)               = (α * r.num) / r.den
Base.:*(r::RationalPoly, α)               = (r.num * α) / r.den

Base.zero(r::RationalPoly) = zero(typeof(r))
Base.zero(::Type{RationalPoly{NT, DT}}) where {NT, DT} = zero(NT) / one(DT)
Base.one(r::RationalPoly) = one(typeof(r))
Base.one(::Type{RationalPoly{NT, DT}}) where {NT, DT} = one(NT) / one(DT)

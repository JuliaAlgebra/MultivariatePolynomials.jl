module MultivariatePolynomialsMooncakeExt
import Mooncake as MC
import ChainRulesCore as CRC
using MultivariatePolynomials
import MultivariatePolynomials as MP
import StarAlgebras as SA
import MutableArithmetics as MA
using Mooncake: Random
const Rules = Base.get_extension(MP, :MultivariatePolynomialsChainRulesCoreExt)
const coefficient_dot = Rules.coefficient_dot
const coefficient_at = Rules.coefficient_at

const CF = Union{Base.IEEEFloat,Complex{<:Base.IEEEFloat}}
const MonomialAlgebra = SA.StarAlgebra{<:MP.Variables{MP.Monomial}}
const FloatPolynomial = SA.AlgebraElement{<:CF,<:MonomialAlgebra}

# A zero primal can have a nonzero coefficient gradient. Its cotangent must
# therefore have its own sparse support, independent of the primal's support.
mutable struct PolynomialTangent{P}
    polynomial::P
end
MC.@foldable function MC.tangent_type(::Type{P}) where {P<:FloatPolynomial}
    return PolynomialTangent{SA.similar_type(P, coefficient_type(P))}
end
MC.@foldable MC.fdata_type(::Type{T}) where {T<:PolynomialTangent} = T
MC.@foldable MC.rdata_type(::Type{<:PolynomialTangent}) = MC.NoRData
MC.@foldable MC.tangent_type(
    ::Type{T},
    ::Type{MC.NoRData},
) where {T<:PolynomialTangent} = T
MC.fdata(t::PolynomialTangent) = t
MC.rdata(::PolynomialTangent) = MC.NoRData()
MC.tangent(t::PolynomialTangent, ::MC.NoRData) = t
function MC.zero_tangent_internal(p::FloatPolynomial, cache::MC.MaybeCache)
    haskey(cache, p) && return cache[p]::MC.tangent_type(typeof(p))
    t = PolynomialTangent(zero(p))
    cache[p] = t
    return t
end
function MC.randn_tangent_internal(
    rng::Random.AbstractRNG,
    p::FloatPolynomial,
    cache::MC.MaybeCache,
)
    haskey(cache, p) && return cache[p]::MC.tangent_type(typeof(p))
    t = PolynomialTangent(
        map_coefficients(_ -> randn(rng, coefficient_type(p)), p),
    )
    cache[p] = t
    return t
end
function MC.increment_internal!!(
    cache::MC.IncCache,
    t::T,
    s::T,
) where {T<:PolynomialTangent}
    (t === s || haskey(cache, t)) && return t
    cache[t] = true
    MA.operate!(+, t.polynomial, s.polynomial)
    return t
end
function MC.set_to_zero_internal!!(
    cache::MC.SetToZeroCache,
    t::PolynomialTangent,
)
    MC._already_tracked!(cache, t) && return t
    MA.operate!(zero, t.polynomial)
    return t
end
function MC._scale_internal(
    cache::MC.MaybeCache,
    a::Float64,
    t::T,
) where {T<:PolynomialTangent}
    haskey(cache, t) && return cache[t]::T
    s = T(coefficient_type(t.polynomial)(a) * t.polynomial)
    cache[t] = s
    return s
end
function MC._dot_internal(
    cache::MC.MaybeCache,
    t::T,
    s::T,
) where {T<:PolynomialTangent}
    key = (t, s)
    haskey(cache, key) && return 0.0
    cache[key] = nothing
    return Float64(real(coefficient_dot(t.polynomial, s.polynomial)))
end
function MC._add_to_primal_internal(
    cache::MC.MaybeCache,
    p::P,
    t::PolynomialTangent,
    unsafe::Bool,
) where {P<:FloatPolynomial}
    key = (p, t, unsafe)
    haskey(cache, key) && return cache[key]::P
    q = convert(P, p + t.polynomial)
    cache[key] = q
    return q
end
function MC.tangent_to_primal_internal!!(
    p::FloatPolynomial,
    t::PolynomialTangent,
    cache::MC.MaybeCache,
)
    haskey(cache, p) && return cache[p]::typeof(p)
    cache[p] = p
    SA.map_coefficients_to!(p, identity, t.polynomial)
    return p
end
function MC.primal_to_tangent_internal!!(
    t::PolynomialTangent,
    p::FloatPolynomial,
    cache::MC.MaybeCache,
)
    haskey(cache, p) && return cache[p]::typeof(t)
    cache[p] = t
    SA.map_coefficients_to!(t.polynomial, identity, p)
    return t
end
function MC.__verify_fdata_value(
    ::IdDict{Any,Nothing},
    ::FloatPolynomial,
    ::PolynomialTangent,
)
    return nothing
end
MC.uninit_fdata(p::FloatPolynomial) = MC.zero_tangent(p)
function MC.TestUtils.populate_address_map_internal(
    m::MC.TestUtils.AddressMap,
    p::FloatPolynomial,
    t::PolynomialTangent,
)
    k, v = pointer_from_objref(p), pointer_from_objref(t)
    haskey(m, k) && (@assert m[k] == v)
    m[k] = v
    return m
end
function MC.TestUtils.has_equal_data_internal(
    t::T,
    s::T,
    ::Bool,
    ::IdDict{Any,Bool},
) where {T<:PolynomialTangent}
    return t.polynomial == s.polynomial
end

# Pullbacks use coefficient polynomials, while a term's only differentiable
# field is its coefficient and scalar arguments have ordinary scalar tangents.
function _cotangent(y::MC.CoDual{<:FloatPolynomial}, ::MC.NoRData)
    return MC.tangent(y).polynomial
end
function _cotangent(y::MC.CoDual{<:AbstractTerm}, dy)
    t = MC.tangent(MC.tangent(y), dy)
    return MC.val(t.fields.coefficient) * monomial(MC.primal(y))
end
_cotangent(y::MC.CoDual{<:Number}, dy) = dy
_scalar_gradient(x::CF, d) = convert(typeof(x), x isa Real ? real(d) : d)
_scalar_gradient(x::Number, d) = MC.NoRData()
function _accumulate!(x::MC.CoDual{<:FloatPolynomial}, d)
    p = MC.primal(x)
    dp = MC.tangent(x).polynomial
    d = d isa Number ? d * one(p) : d
    for t in terms(d)
        m = monomial(t)
        if all(v -> iszero(degree(m, v)) || v in variables(p), variables(m))
            c = _scalar_gradient(zero(coefficient_type(p)), coefficient(t))
            # Translate the key to the input's fixed parent after dropping
            # variables whose exponents are zero.
            k = map(Base.Fix1(degree, m), variables(p))
            MA.operate!(+, dp, SA.Term(parent(p), k, c))
        end
    end
    return MC.NoRData()
end
function _accumulate!(x::MC.CoDual{<:AbstractTerm}, d)
    c = coefficient(MC.primal(x))
    dc = _scalar_gradient(c, coefficient_at(d, monomial(MC.primal(x))))
    return MC.increment_field!!(
        MC.zero_rdata(MC.primal(x)),
        dc,
        Val(:coefficient),
    )
end
_accumulate!(x::MC.CoDual{<:Number}, d) = _scalar_gradient(MC.primal(x), d)
function _accumulate!(x::MC.CoDual{<:AbstractVector{<:CF}}, d)
    px, dx = MC.arrayify(x)
    dx .+= _scalar_gradient.(px, d)
    return MC.NoRData()
end
function _accumulate(x, d)
    if d isa CRC.AbstractZero || (
        MC.primal(x) isa MP._APL &&
        MC.tangent_type(coefficient_type(MC.primal(x))) === MC.NoTangent
    )
        return MC.zero_rdata(MC.primal(x))
    end
    return _accumulate!(x, d)
end
function _accumulate!(x::MC.CoDual{<:AbstractMonomialLike}, d)
    return MC.zero_rdata(MC.primal(x))
end

function _rrule(f, args...)
    y, pb = CRC.rrule(MC.primal(f), map(MC.primal, args)...)
    out = MC.zero_fcodual(y)
    if MC.tangent_type(typeof(y)) === MC.NoTangent ||
       (y isa MP._APL && MC.tangent_type(coefficient_type(y)) === MC.NoTangent)
        return out, MC.NoPullback(f, args...)
    end
    if length(args) == 1 && y isa FloatPolynomial && y === MC.primal(only(args))
        return only(args), MC.NoPullback(f, args...)
    end
    function pullback(dy)
        ds = Base.tail(pb(_cotangent(out, dy)))
        return MC.NoRData(), map(_accumulate, args, ds)...
    end
    return out, pullback
end
for op in (+, -, *)
    for types in ((MP._APL, MP._APL), (MP._APL, Number), (Number, MP._APL))
        A, B = types
        @eval begin
            MC.@is_primitive MC.DefaultCtx MC.ReverseMode Tuple{
                typeof($op),
                $A,
                $B,
            }
            MC.rrule!!(
                f::MC.CoDual{typeof($op)},
                a::MC.CoDual{<:$A},
                b::MC.CoDual{<:$B},
            ) = _rrule(f, a, b)
        end
    end
end
for op in (+, -, polynomial)
    @eval begin
        MC.@is_primitive MC.DefaultCtx MC.ReverseMode Tuple{typeof($op),MP._APL}
        MC.rrule!!(f::MC.CoDual{typeof($op)}, a::MC.CoDual{<:MP._APL}) =
            _rrule(f, a)
    end
end
for (op, A, B) in (
    (/, MP._APL, Number),
    (^, MP._APL, Integer),
    (differentiate, MP._APL, AbstractVariable),
    (coefficient, MP._APL, AbstractMonomialLike),
    (term, Number, AbstractMonomialLike),
    (
        polynomial,
        AbstractVector{<:Number},
        AbstractVector{<:AbstractMonomialLike},
    ),
)
    @eval begin
        MC.@is_primitive MC.DefaultCtx MC.ReverseMode Tuple{typeof($op),$A,$B}
        MC.rrule!!(
            f::MC.CoDual{typeof($op)},
            a::MC.CoDual{<:$A},
            b::MC.CoDual{<:$B},
        ) = _rrule(f, a, b)
    end
end
MC.@is_primitive MC.DefaultCtx MC.ReverseMode Tuple{
    typeof(coefficient),
    AbstractTerm,
}
function MC.rrule!!(
    f::MC.CoDual{typeof(coefficient)},
    a::MC.CoDual{<:AbstractTerm},
)
    return _rrule(f, a)
end
end

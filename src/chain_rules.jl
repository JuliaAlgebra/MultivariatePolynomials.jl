# The publlback depends on the scalar product on the polynomials
# With the scalar product `LinearAlgebra.dot(p, q) = p * q`, there is no pullback for `differentiate`
# With the scalar product `_dot(p, q)` of `test/chain_rules.jl`, there is a pullback for `differentiate`
# and the pullback for `*` changes.
# We give the one for the scalar product `_dot`.

import ChainRulesCore

ChainRulesCore.@scalar_rule +(x::APL) true
ChainRulesCore.@scalar_rule -(x::APL) -1

ChainRulesCore.@scalar_rule +(x::APL, y::APL) (true, true)
function plusconstant1_pullback(Δ)
    return ChainRulesCore.NoTangent(), Δ, coefficient(Δ, constantmonomial(Δ))
end
function ChainRulesCore.rrule(::typeof(plusconstant), p::APL, α)
    return plusconstant(p, α), plusconstant1_pullback
end
function plusconstant2_pullback(Δ)
    return ChainRulesCore.NoTangent(), coefficient(Δ, constantmonomial(Δ)), Δ
end
function ChainRulesCore.rrule(::typeof(plusconstant), α, p::APL)
    return plusconstant(α, p), plusconstant2_pullback
end
ChainRulesCore.@scalar_rule -(x::APL, y::APL) (true, -1)

function ChainRulesCore.frule((_, Δp, Δq), ::typeof(*), p::APL, q::APL)
    return p * q, MA.add_mul!!(p * Δq, q, Δp)
end

function _mult_pullback(op::F, ts, p, Δ) where {F<:Function}
    for t in terms(p)
        c = coefficient(t)
        m = monomial(t)
        for δ in Δ
            if divides(m, δ)
                coef = op(c, coefficient(δ))
                mono = _div(monomial(δ), m)
                push!(ts, term(coef, mono))
            end
        end
    end
    return polynomial(ts)
end
function adjoint_mult_left(p, Δ)
    ts = MA.promote_operation(*, MA.promote_operation(adjoint, termtype(p)), termtype(Δ))[]
    return _mult_pullback(ts, p, Δ) do c, d
        c' * d
    end
end
function adjoint_mult_right(p, Δ)
    ts = MA.promote_operation(*, termtype(Δ), MA.promote_operation(adjoint, termtype(p)))[]
    return _mult_pullback(ts, p, Δ) do c, d
        d * c'
    end
end

function ChainRulesCore.rrule(::typeof(*), p::APL, q::APL)
    function times_pullback2(ΔΩ̇)
        # This is for the scalar product `_dot`:
        return (ChainRulesCore.NoTangent(), adjoint_mult_right(q, ΔΩ̇), adjoint_mult_left(p, ΔΩ̇))
        # For the scalar product `dot`, it would be instead:
        return (ChainRulesCore.NoTangent(), ΔΩ̇ * q', p' * ΔΩ̇)
    end
    return p * q, times_pullback2
end

function ChainRulesCore.rrule(::typeof(multconstant), α, p::APL)
    function times_pullback2(ΔΩ̇)
        # TODO we could make it faster, don't need to compute `Δα` entirely if we only care about the constant term.
        Δα = adjoint_mult_right(p, ΔΩ̇)
        return (ChainRulesCore.NoTangent(), coefficient(Δα, constantmonomial(Δα)), α' * ΔΩ̇)
    end
    return multconstant(α, p), times_pullback2
end

function ChainRulesCore.rrule(::typeof(multconstant), p::APL, α)
    function times_pullback2(ΔΩ̇)
        # TODO we could make it faster, don't need to compute `Δα` entirely if we only care about the constant term.
        Δα = adjoint_mult_left(p, ΔΩ̇)
        return (ChainRulesCore.NoTangent(), ΔΩ̇ * α', coefficient(Δα, constantmonomial(Δα)))
    end
    return multconstant(p, α), times_pullback2
end

notangent3(Δ) = ChainRulesCore.NoTangent(), ChainRulesCore.NoTangent(), ChainRulesCore.NoTangent()
function ChainRulesCore.rrule(::typeof(^), mono::AbstractMonomialLike, i::Integer)
    return mono^i, notangent3
end

function ChainRulesCore.frule((_, Δp, _), ::typeof(differentiate), p, x)
    return differentiate(p, x), differentiate(Δp, x)
end
# This is for the scalar product `_dot`, there is no pullback for the scalar product `dot`
function differentiate_pullback(Δdpdx, x)
    return ChainRulesCore.NoTangent(), x * differentiate(x * Δdpdx, x), ChainRulesCore.NoTangent()
end
function ChainRulesCore.rrule(::typeof(differentiate), p, x)
    dpdx = differentiate(p, x)
    return dpdx, Base.Fix2(differentiate_pullback, x)
end

function coefficient_pullback(Δ, m::AbstractMonomialLike)
    return ChainRulesCore.NoTangent(), polynomial(term(Δ, m)), ChainRulesCore.NoTangent()
end
function ChainRulesCore.rrule(::typeof(coefficient), p::APL, m::AbstractMonomialLike)
    return coefficient(p, m), Base.Fix2(coefficient_pullback, m)
end

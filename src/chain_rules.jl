import ChainRulesCore

ChainRulesCore.@scalar_rule +(x::APL, y::APL) (true, true)
ChainRulesCore.@scalar_rule -(x::APL, y::APL) (true, -1)

function ChainRulesCore.frule((_, Δp, Δq), ::typeof(*), p::APL, q::APL)
    return p * q, MA.add_mul!!(p * Δq, q, Δp)
end
function ChainRulesCore.rrule(::typeof(*), p::APL, q::APL)
    function times_pullback2(ΔΩ̇)
        #ΔΩ = ChainRulesCore.unthunk(Ω̇)
        #return (ChainRulesCore.NoTangent(), ChainRulesCore.ProjectTo(p)(ΔΩ * q'), ChainRulesCore.ProjectTo(q)(p' * ΔΩ))
        return (ChainRulesCore.NoTangent(), ΔΩ̇ * q', p' * ΔΩ̇)
    end
    return p * q, times_pullback2
end

function ChainRulesCore.frule((_, Δp, _), ::typeof(differentiate), p, x)
    return differentiate(p, x), differentiate(Δp, x)
end
function pullback_differentiate_polynomial(Δdpdx, x)
    return ChainRulesCore.NoTangent(), x * differentiate(x * Δdpdx, x), ChainRulesCore.NoTangent()
end
function ChainRulesCore.rrule(::typeof(differentiate), p::APL, x)
    dpdx = differentiate(p, x)
    return dpdx, Base.Fix2(pullback_differentiate_polynomial, x)
end

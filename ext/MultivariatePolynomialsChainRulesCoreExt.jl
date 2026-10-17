module MultivariatePolynomialsChainRulesCoreExt

import ChainRulesCore
using MultivariatePolynomials
using MultivariatePolynomials: _APL, MA

ChainRulesCore.@scalar_rule +(x::_APL) true
ChainRulesCore.@scalar_rule -(x::_APL) -1

ChainRulesCore.@scalar_rule +(x::_APL, y::_APL) (true, true)
ChainRulesCore.@scalar_rule -(x::_APL, y::_APL) (true, -1)

function ChainRulesCore.frule((_, Δp, Δq), ::typeof(*), p::_APL, q::_APL)
    return p * q, MA.add_mul!!(p * Δq, q, Δp)
end
function ChainRulesCore.rrule(::typeof(*), p::_APL, q::_APL)
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
function pullback(Δdpdx, x)
    return ChainRulesCore.NoTangent(),
    x * differentiate(x * Δdpdx, x),
    ChainRulesCore.NoTangent()
end
function ChainRulesCore.rrule(::typeof(differentiate), p, x)
    dpdx = differentiate(p, x)
    return dpdx, Base.Fix2(pullback, x)
end

end

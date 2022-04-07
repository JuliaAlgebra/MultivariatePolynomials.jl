import ChainRulesCore

ChainRulesCore.@scalar_rule +(x::APL, y::APL) (true, true)
function ChainRulesCore.frule((_, Δp, _), ::typeof(differentiate), p, x)
    return differentiate(p, x), differentiate(Δp, x)
end
function ChainRulesCore.rrule(::typeof(differentiate), p, x)
    dpdx = differentiate(p, x)
    function pullback(Δdpdx)
        return ChainRulesCore.NoTangent(), x * differentiate(x * Δdpdx, x), ChainRulesCore.NoTangent()
    end
    return dpdx, pullback
end

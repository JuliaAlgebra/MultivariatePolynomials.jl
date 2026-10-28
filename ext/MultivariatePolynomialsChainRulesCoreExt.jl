module MultivariatePolynomialsChainRulesCoreExt

import ChainRulesCore
import ChainRulesCore as CRC
import MultivariatePolynomials as MP
using MultivariatePolynomials
using MultivariatePolynomials: _APL, MA

ChainRulesCore.@scalar_rule +(x::_APL) true
ChainRulesCore.@scalar_rule -(x::_APL) -1

ChainRulesCore.@scalar_rule +(x::_APL, y::_APL) (true, true)
ChainRulesCore.@scalar_rule -(x::_APL, y::_APL) (true, -1)

function ChainRulesCore.frule((_, Δp, Δq), ::typeof(*), p::_APL, q::_APL)
    return p * q, MA.add_mul!!(p * Δq, q, Δp)
end

function ChainRulesCore.frule((_, Δp, _), ::typeof(differentiate), p, x)
    return differentiate(p, x), differentiate(Δp, x)
end
include("coefficient_rules.jl")

# Share the coefficient-space pullbacks with the Mooncake extension.
for op in (+, -, *)
    for (A, B) in ((_APL, Number), (Number, _APL))
        @eval CRC.rrule(::typeof($op), p::$A, q::$B) = _rrule($op, p, q)
    end
end
for (op, A, B) in (
    (*, _APL, _APL),
    (/, _APL, Number),
    (^, _APL, Integer),
    (differentiate, _APL, AbstractVariable),
    (coefficient, _APL, AbstractMonomialLike),
    (term, Number, AbstractMonomialLike),
    (
        polynomial,
        AbstractVector{<:Number},
        AbstractVector{<:AbstractMonomialLike},
    ),
)
    @eval CRC.rrule(::typeof($op), p::$A, q::$B) = _rrule($op, p, q)
end
CRC.rrule(::typeof(coefficient), p::AbstractTerm) = _rrule(coefficient, p)
CRC.rrule(::typeof(polynomial), p::_APL) = _rrule(polynomial, p)

end

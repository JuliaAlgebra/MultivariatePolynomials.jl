# Fast univariate `gcd` used as basis for multivariate `gcd`.
# Note that we cannot use `Polynomials.jl` as it assumes that the coefficients
# are `<:Number` while multivariate `gcd` creates univariate polynomials
# for which the coefficients are polynomials in the rest of the variables
struct UnivariateMonomial <: AbstractMonomial
    exponent::Int
end

const UnivariateTerm{T} = Term{T, UnivariateMonomial}
const UnivariatePolynomial{T} = Polynomial{T, UnivariateMonomial}

function MA.mutable_operate!(::typeof(univariate_rem), f::APL, g::APL)
    ltg = leadingterm(g)
    ltf = leadingterm(f)
    while divides(monomial(ltg), ltf)
        MA.mutable_operate!(multconstant, f, coefficient(ltg))
    end
end
function univariate_gcd!!(p1, p2)
end

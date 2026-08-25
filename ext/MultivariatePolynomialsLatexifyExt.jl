module MultivariatePolynomialsLatexifyExt
using MultivariatePolynomials
using Latexify
@latexrecipe function f(m::AbstractMonomial)
    return Expr(:call, :*, (:($(Symbol(effvar)) ^ $e) for (effvar, e) = zip(effective_variables(m), exponents(m)))...)
end
end
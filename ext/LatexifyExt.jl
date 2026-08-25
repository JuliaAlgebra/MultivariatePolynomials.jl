module LatexifyExt
using MultivariatePolynomials
using Latexify
@latexrecipe function f(m::AbstractMonomialLike)
    return Expr(:call, :*, (:($(Symbol(effvar)) ^ $e) for (effvar, e) = zip(effective_variables(m), exponents(m)))...)
end

@latexrecipe function f(t::AbstractTermLike)
    return Expr(:call, :*, coefficient(t), monomial(t))
end

@latexrecipe function f(p::AbstractPolynomialLike)
    return Expr(:call, :+, terms(p)...)
end
end
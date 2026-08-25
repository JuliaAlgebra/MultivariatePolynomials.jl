module MultivariatePolynomialsLatexifyExt
using MultivariatePolynomials
using Latexify
@latexrecipe function f(m::AbstractMonomial)
    return Expr(:call, :*, (:($(Symbol(effvar)) ^ $e) for (effvar, e) = zip(effective_variables(m), exponents(m)))...)
end

@latexrecipe function f(t::AbstractTerm)
           return Expr(:call, :*, coefficient(t), monomial(t))
       end

@latexrecipe function f(p::Polynomial)
    return Expr(:call, :+, terms(p)...)
end
end
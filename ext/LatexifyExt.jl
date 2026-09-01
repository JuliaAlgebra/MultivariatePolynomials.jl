module LatexifyExt
using MultivariatePolynomials
using Latexify
@latexrecipe function f(m::AbstractMonomialLike)
    operation := :*
    mult_symbol --> ""
    vars = variables(m)
    exps = exponents(m)
    factors = ((e == 1) ? Symbol(var) : :($(Symbol(var)) ^ $e) for (var, e) = zip(vars, exps) if e != 0)
    return Expr(:call, :*, factors...)
end

@latexrecipe function f(t::AbstractTermLike)
    return Expr(:call, :*, coefficient(t), monomial(t))
end

@latexrecipe function f(p::AbstractPolynomialLike)
    return Expr(:call, :+, terms(p)...)
end
end
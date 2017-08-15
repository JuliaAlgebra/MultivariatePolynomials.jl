var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#MultivariatePolynomials-1",
    "page": "Introduction",
    "title": "MultivariatePolynomials",
    "category": "section",
    "text": "MultivariatePolynomials.jl is an implementation independent library for manipulating multivariate polynomials. It defines abstract types and an API for multivariate monomials, terms, polynomials, moments and measures and gives default implementation for common operations on them using the API. If you want to manipulate multivariate polynomials easily and efficiently while being able to easily switch between different implementations, this library is exactly what you are looking for.Supported operations are : basic arithmetic, rational polynomials, differentiation and evaluation/substitution, division and duality operations between polynomials and moments. There is also support for solving systems of equations (soon!) and building (semi)algebraic sets.Currently, the following implementations are available:TypedPolynomials\nDynamicPolynomialsPages = [\"apireference.md\"]\nDepth = 3"
},

{
    "location": "apireference.html#",
    "page": "Reference",
    "title": "Reference",
    "category": "page",
    "text": "CurrentModule = MultivariatePolynomials"
},

{
    "location": "apireference.html#API-1",
    "page": "Reference",
    "title": "API",
    "category": "section",
    "text": ""
},

{
    "location": "apireference.html#MultivariatePolynomials.name",
    "page": "Reference",
    "title": "MultivariatePolynomials.name",
    "category": "Function",
    "text": "name(v::AbstractVariable)::AbstractString\n\nReturns the name of a variable.\n\n\n\n"
},

{
    "location": "apireference.html#Variables-1",
    "page": "Reference",
    "title": "Variables",
    "category": "section",
    "text": "name"
},

{
    "location": "apireference.html#MultivariatePolynomials.nvariables",
    "page": "Reference",
    "title": "MultivariatePolynomials.nvariables",
    "category": "Function",
    "text": "nvariables(p::AbstractPolynomialLike)\n\nReturns the number of variables in p, i.e. length(variables(p)). It could be more than the number of variables with nonzero exponent (see the Examples section of variables).\n\nExamples\n\nCalling nvariables(x^2*y) should return at least 2 and calling nvariables(x) should return at least 1.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.variables",
    "page": "Reference",
    "title": "MultivariatePolynomials.variables",
    "category": "Function",
    "text": "variables(p::AbstractPolynomialLike)\n\nReturns the tuple of the variables of p in decreasing order. It could contain variables of zero degree, see the example section.\n\nExamples\n\nCalling variables(x^2*y) should return (x, y) and calling variables(x) should return (x,). Note that the variables of m does not necessarily have nonzero exponent. For instance, variables([x^2*y, y*z][1]) is usually (x, y, z) since the two monomials have been promoted to a common type.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.exponents",
    "page": "Reference",
    "title": "MultivariatePolynomials.exponents",
    "category": "Function",
    "text": "exponents(t::AbstractTermLike)\n\nReturns the exponent of the variables in the monomial of the term t.\n\nExamples\n\nCalling exponents(x^2*y) should return (2, 1).\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.exponent",
    "page": "Reference",
    "title": "MultivariatePolynomials.exponent",
    "category": "Function",
    "text": "exponent(t::AbstractTermLike, var::AbstractVariable)\n\nReturns the exponent of the variable var in the monomial of the term t.\n\nExamples\n\nCalling exponent(x^2*y, x) should return 2 and calling exponent(x^2*y, y) should return 1.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.deg",
    "page": "Reference",
    "title": "MultivariatePolynomials.deg",
    "category": "Function",
    "text": "deg(t::AbstractTermLike)\n\nReturns the total degree of the monomial of the term t, i.e. sum(exponents(t)).\n\nExamples\n\nCalling deg(x^2*y) should return 3 which is 2 + 1.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.isconstant",
    "page": "Reference",
    "title": "MultivariatePolynomials.isconstant",
    "category": "Function",
    "text": "isconstant(m::AbstractMonomialLike)\n\nReturns whether the monomial is constant.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.constantmonomial",
    "page": "Reference",
    "title": "MultivariatePolynomials.constantmonomial",
    "category": "Function",
    "text": "constantmonomial(p::AbstractPolynomialType)\n\nReturns a constant monomial of the monomial type of p with the same variables as p.\n\nconstantmonomial(::Type{PT}) where {PT<:AbstractPolynomialType}\n\nReturns a constant monomial of the monomial type of a polynomial of type PT.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.divides",
    "page": "Reference",
    "title": "MultivariatePolynomials.divides",
    "category": "Function",
    "text": "divides(t1::AbstractTermLike, t2::AbstractTermLike)\n\nReturns whether monomial(t1) divides monomial(t2).\n\nExamples\n\nCalling divides(2x^2y, 3xy) should return false because x^2y does not divide xy since x has a degree 2 in x^2y which is greater than the degree of x on xy. However, calling divides(3xy, 2x^2y) should return true.\n\n\n\n"
},

{
    "location": "apireference.html#Monomials-1",
    "page": "Reference",
    "title": "Monomials",
    "category": "section",
    "text": "nvariables\nvariables\nexponents\nexponent\ndeg\nisconstant\nconstantmonomial\ndivides"
},

{
    "location": "apireference.html#MultivariatePolynomials.term",
    "page": "Reference",
    "title": "MultivariatePolynomials.term",
    "category": "Function",
    "text": "term(p::AbstractPolynomialLike)\n\nConverts the polynomial p to a term. When applied on a polynomial, it throws an error if it has more than one term. When applied to a term, it is the identity and does not copy it. When applied to a monomial, it create a term of type AbstractTerm{Int}.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.coefficient",
    "page": "Reference",
    "title": "MultivariatePolynomials.coefficient",
    "category": "Function",
    "text": "coefficient(t::AbstractTermLike)\n\nReturns the coefficient of the term t.\n\nExamples\n\nCalling coefficient on 4x^2y should return 4.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.monomial",
    "page": "Reference",
    "title": "MultivariatePolynomials.monomial",
    "category": "Function",
    "text": "monomial(t::AbstractTermLike)\n\nReturns the monomial of the term t.\n\nExamples\n\nCalling monomial on 4x^2y should return x^2y.\n\n\n\n"
},

{
    "location": "apireference.html#Terms-1",
    "page": "Reference",
    "title": "Terms",
    "category": "section",
    "text": "term\ncoefficient\nmonomial"
},

{
    "location": "apireference.html#MultivariatePolynomials.terms",
    "page": "Reference",
    "title": "MultivariatePolynomials.terms",
    "category": "Function",
    "text": "terms(p::AbstractPolynomialLike)\n\nReturns an iterator over the nonzero terms of the polynomial p sorted in the decreasing monomial order.\n\nExamples\n\nCalling terms on 4x^2y + xy + 2x should return an iterator of 4x^2y xy 2x.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.monomials",
    "page": "Reference",
    "title": "MultivariatePolynomials.monomials",
    "category": "Function",
    "text": "monomials(p::AbstractPolynomialLike)\n\nReturns an iterator over the monomials of p of the nonzero terms of the polynomial sorted in the decreasing order.\n\nmonomials(vars::Tuple, degs::AbstractVector{Int}, filter::Function = m -> true)\n\nBuilds the vector of all the monovec m with variables vars such that the degree deg(m) is in degs and filter(m) is true.\n\nExamples\n\nCalling monomials on 4x^2y + xy + 2x should return an iterator of x^2y xy x.\n\nCalling monomials((x, y), [1, 3], m -> exponent(m, y) != 1) should return [x^3, x*y^2, y^3, x] where x^2*y and y have been excluded by the filter.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.mindeg",
    "page": "Reference",
    "title": "MultivariatePolynomials.mindeg",
    "category": "Function",
    "text": "mindeg(p::AbstractPolynomialLike)\n\nReturns the minimal total degree of the monomials of p, i.e. minimum(deg, terms(p)).\n\nExamples\n\nCalling mindeg on on 4x^2y + xy + 2x should return 1.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.maxdeg",
    "page": "Reference",
    "title": "MultivariatePolynomials.maxdeg",
    "category": "Function",
    "text": "maxdeg(p::AbstractPolynomialLike)\n\nReturns the maximal total degree of the monomials of p, i.e. maximum(deg, terms(p)).\n\nExamples\n\nCalling maxdeg on on 4x^2y + xy + 2x should return 3.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.extdeg",
    "page": "Reference",
    "title": "MultivariatePolynomials.extdeg",
    "category": "Function",
    "text": "extdeg(p::AbstractPolynomialLike)\n\nReturns the extremal total degrees of the monomials of p, i.e. (mindeg(p), maxdeg(p)).\n\nExamples\n\nCalling extdeg on on 4x^2y + xy + 2x should return (1, 3).\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.leadingterm",
    "page": "Reference",
    "title": "MultivariatePolynomials.leadingterm",
    "category": "Function",
    "text": "leadingterm(p::AbstractPolynomialLike)\n\nReturns the coefficient of the leading term, i.e. first(terms(p)).\n\nExamples\n\nCalling leadingterm on 4x^2y + xy + 2x should return 4x^2y.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.leadingcoefficient",
    "page": "Reference",
    "title": "MultivariatePolynomials.leadingcoefficient",
    "category": "Function",
    "text": "leadingcoefficient(p::AbstractPolynomialLike)\n\nReturns the coefficient of the leading term of p, i.e. coefficient(leadingterm(p)).\n\nExamples\n\nCalling leadingcoefficient on 4x^2y + xy + 2x should return 4 and calling it on 0 should return 0.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.leadingmonomial",
    "page": "Reference",
    "title": "MultivariatePolynomials.leadingmonomial",
    "category": "Function",
    "text": "leadingmonomial(p::AbstractPolynomialLike)\n\nReturns the monomial of the leading term of p, i.e. monomial(leadingterm(p)) or first(monomials(p)).\n\nExamples\n\nCalling leadingmonomial on 4x^2y + xy + 2x should return x^2y.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.removeleadingterm",
    "page": "Reference",
    "title": "MultivariatePolynomials.removeleadingterm",
    "category": "Function",
    "text": "removeleadingterm(p::AbstractPolynomialLike)\n\nReturns a polynomial with the leading term removed in the polynomial p.\n\nExamples\n\nCalling removeleadingterm on 4x^2y + xy + 2x should return xy + 2x.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.removemonomials",
    "page": "Reference",
    "title": "MultivariatePolynomials.removemonomials",
    "category": "Function",
    "text": "Returns a polynomial with the terms having their monomial in the monomial vector mv removed in the polynomial p.\n\nExamples\n\nCalling removemonomials(4x^2*y + x*y + 2x, [x*y]) should return 4x^2*y + 2x.\n\n\n\n"
},

{
    "location": "apireference.html#Polynomials-1",
    "page": "Reference",
    "title": "Polynomials",
    "category": "section",
    "text": "terms\nmonomials\nmindeg\nmaxdeg\nextdeg\nleadingterm\nleadingcoefficient\nleadingmonomial\nremoveleadingterm\nremovemonomials"
},

{
    "location": "apireference.html#MultivariatePolynomials.monovec",
    "page": "Reference",
    "title": "MultivariatePolynomials.monovec",
    "category": "Function",
    "text": "monovec(X::AbstractVector{MT}) where {MT<:AbstractMonomialLike}\n\nReturns the vector of monomials X in decreasing order and without any duplicates.\n\nExamples\n\nCalling monovec on xy x xy x^2y x should return x^2y xy x.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.monovectype",
    "page": "Reference",
    "title": "MultivariatePolynomials.monovectype",
    "category": "Function",
    "text": "monovectype(X::AbstractVector{MT}) where {MT<:AbstractMonomialLike}\n\nReturns the return type of monovec.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.sortmonovec",
    "page": "Reference",
    "title": "MultivariatePolynomials.sortmonovec",
    "category": "Function",
    "text": "sortmonovec(X::AbstractVector{MT}) where {MT<:AbstractMonomialLike}\n\nReturns σ, the orders in which one must take the monomials in X to make them sorted and without any duplicate and the sorted vector of monomials, i.e. it returns (σ, X[σ]).\n\nExamples\n\nCalling sortmonovec on xy x xy x^2y x should return (4 1 2 x^2y xy x).\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.mergemonovec",
    "page": "Reference",
    "title": "MultivariatePolynomials.mergemonovec",
    "category": "Function",
    "text": "mergemonovec{MT<:AbstractMonomialLike, MVT<:AbstractVector{MT}}(X::AbstractVector{MVT}}\n\nReturns the vector of monomials in the entries of X in decreasing order and without any duplicates, i.e. monovec(vcat(X...))\n\nExamples\n\nCalling mergemonovec on xy x xy x^2y x should return x^2y xy x.\n\n\n\n"
},

{
    "location": "apireference.html#Monomial-Vectors-1",
    "page": "Reference",
    "title": "Monomial Vectors",
    "category": "section",
    "text": "monovec\nmonovectype\nsortmonovec\nmergemonovec"
},

]}

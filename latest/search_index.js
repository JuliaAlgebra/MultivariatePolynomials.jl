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
    "location": "apireference.html#Terms-1",
    "page": "Reference",
    "title": "Terms",
    "category": "section",
    "text": "term\nzeroterm\ncoefficient\nmonomial\nexponent\ndeg\nisconstant\ndivides"
},

{
    "location": "apireference.html#Polynomials-1",
    "page": "Reference",
    "title": "Polynomials",
    "category": "section",
    "text": "terms\nmonomials\nmindeg\nmaxdeg\nextdeg\nleadingterm\nleadingcoefficient\nleadingmonomial\nremoveleadingterm\nremovemonomials\nvars\nnvars"
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

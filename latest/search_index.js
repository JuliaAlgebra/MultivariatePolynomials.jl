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
    "location": "apireference.html#MultivariatePolynomials.variable",
    "page": "Reference",
    "title": "MultivariatePolynomials.variable",
    "category": "Function",
    "text": "variable(p::AbstractPolynomialLike)\n\nConverts p to a variable. Throws an error if it is not possible.\n\nExamples\n\nCalling variable(x^2 + x - x^2) should return the variable x and calling variable(1.0y) should return the variable y however calling variable(2x) or variable(x + y) should throw an error.\n\nNote\n\nThis operation is not type stable for the TypedPolynomials implementation if nvariables(p) > 1 but is type stable for DynamicPolynomials.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.name",
    "page": "Reference",
    "title": "MultivariatePolynomials.name",
    "category": "Function",
    "text": "name(v::AbstractVariable)::AbstractString\n\nReturns the name of a variable.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.similarvariable",
    "page": "Reference",
    "title": "MultivariatePolynomials.similarvariable",
    "category": "Function",
    "text": "similarvariable(p::AbstractPolynomialLike, variable::Type{Val{V}})\n\nCreates a new variable V based upon the the given source polynomial.\n\nsimilarvariable(p::AbstractPolynomialLike, v::Symbol)\n\nCreates a new variable based upon the given source polynomial and the given symbol v. Note that this can lead to type instabilities.\n\nExamples\n\nCalling similarvariable(typedpoly, Val{:x}) on a polynomial created with TypedPolynomials results in TypedPolynomials.Variable{:x}.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.@similarvariable",
    "page": "Reference",
    "title": "MultivariatePolynomials.@similarvariable",
    "category": "Macro",
    "text": "@similarvariable(p::AbstractPolynomialLike, variable)\n\nCalls similarvariable(p, Val{variable}) and binds the result to a variable with the same name.\n\nExamples\n\nCalling @similarvariable typedpoly x on a polynomial created with TypedPolynomials binds TypedPolynomials.Variable{:x} to the variable x.\n\n\n\n"
},

{
    "location": "apireference.html#Variables-1",
    "page": "Reference",
    "title": "Variables",
    "category": "section",
    "text": "variable\nname\nsimilarvariable\n@similarvariable"
},

{
    "location": "apireference.html#MultivariatePolynomials.monomialtype",
    "page": "Reference",
    "title": "MultivariatePolynomials.monomialtype",
    "category": "Function",
    "text": "monomialtype(p::AbstractPolynomialLike)\n\nReturn the type of the monomials of p.\n\ntermtype(::Type{PT}) where PT<:AbstractPolynomialLike\n\nReturns the type of the monomials of a polynomial of type PT.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.variables",
    "page": "Reference",
    "title": "MultivariatePolynomials.variables",
    "category": "Function",
    "text": "variables(p::AbstractPolynomialLike)\n\nReturns the tuple of the variables of p in decreasing order. It could contain variables of zero degree, see the example section.\n\nExamples\n\nCalling variables(x^2*y) should return (x, y) and calling variables(x) should return (x,). Note that the variables of m does not necessarily have nonzero degree. For instance, variables([x^2*y, y*z][1]) is usually (x, y, z) since the two monomials have been promoted to a common type.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.nvariables",
    "page": "Reference",
    "title": "MultivariatePolynomials.nvariables",
    "category": "Function",
    "text": "nvariables(p::AbstractPolynomialLike)\n\nReturns the number of variables in p, i.e. length(variables(p)). It could be more than the number of variables with nonzero degree (see the Examples section of variables).\n\nExamples\n\nCalling nvariables(x^2*y) should return at least 2 and calling nvariables(x) should return at least 1.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.exponents",
    "page": "Reference",
    "title": "MultivariatePolynomials.exponents",
    "category": "Function",
    "text": "exponents(t::AbstractTermLike)\n\nReturns the exponent of the variables in the monomial of the term t.\n\nExamples\n\nCalling exponents(x^2*y) should return (2, 1).\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.degree",
    "page": "Reference",
    "title": "MultivariatePolynomials.degree",
    "category": "Function",
    "text": "degree(t::AbstractTermLike)\n\nReturns the total degree of the monomial of the term t, i.e. sum(exponents(t)).\n\ndegree(t::AbstractTermLike, v::AbstractVariable)\n\nReturns the exponent of the variable v in the monomial of the term t.\n\nExamples\n\nCalling degree(x^2*y) should return 3 which is 2 + 1. Calling degree(x^2*y, x) should return 2 and calling degree(x^2*y, y) should return 1.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.isconstant",
    "page": "Reference",
    "title": "MultivariatePolynomials.isconstant",
    "category": "Function",
    "text": "isconstant(t::AbstractTermLike)\n\nReturns whether the monomial of t is constant.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.powers",
    "page": "Reference",
    "title": "MultivariatePolynomials.powers",
    "category": "Function",
    "text": "powers(t::AbstractTermLike)\n\nReturns an tuple of the powers of the monomial of t.\n\nExamples\n\nCalling powers(3x^4*y) should return((x, 4), (y, 1))`.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.divides",
    "page": "Reference",
    "title": "MultivariatePolynomials.divides",
    "category": "Function",
    "text": "divides(t1::AbstractTermLike, t2::AbstractTermLike)\n\nReturns whether monomial(t1) divides monomial(t2).\n\nExamples\n\nCalling divides(2x^2y, 3xy) should return false because x^2y does not divide xy since x has a degree 2 in x^2y which is greater than the degree of x on xy. However, calling divides(3xy, 2x^2y) should return true.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.constantmonomial",
    "page": "Reference",
    "title": "MultivariatePolynomials.constantmonomial",
    "category": "Function",
    "text": "constantmonomial(p::AbstractPolynomialType)\n\nReturns a constant monomial of the monomial type of p with the same variables as p.\n\nconstantmonomial(::Type{PT}) where {PT<:AbstractPolynomialType}\n\nReturns a constant monomial of the monomial type of a polynomial of type PT.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.mapexponents",
    "page": "Reference",
    "title": "MultivariatePolynomials.mapexponents",
    "category": "Function",
    "text": "mapexponents(f, m1::AbstractMonomialLike, m2::AbstractMonomialLike)\n\nIf m_1 = prod x^alpha_i and m_2 = prod x^beta_i then it returns the monomial m = prod x^f(alpha_i beta_i).\n\nExamples\n\nThe multiplication m1 * m2 is equivalent to mapexponents(+, m1, m2), the unsafe division _div(m1, m2) is equivalent to mapexponents(-, m1, m2), gcd(m1, m2) is equivalent to mapexponents(min, m1, m2), lcm(m1, m2) is equivalent to mapexponents(max, m1, m2).\n\n\n\n"
},

{
    "location": "apireference.html#Monomials-1",
    "page": "Reference",
    "title": "Monomials",
    "category": "section",
    "text": "monomialtype\nvariables\nnvariables\nexponents\ndegree\nisconstant\npowers\ndivides\nconstantmonomial\nmapexponents"
},

{
    "location": "apireference.html#MultivariatePolynomials.term",
    "page": "Reference",
    "title": "MultivariatePolynomials.term",
    "category": "Function",
    "text": "term(p::AbstractPolynomialLike)\n\nConverts the polynomial p to a term. When applied on a polynomial, it throws an error if it has more than one term. When applied to a term, it is the identity and does not copy it. When applied to a monomial, it create a term of type AbstractTerm{Int}.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.termtype",
    "page": "Reference",
    "title": "MultivariatePolynomials.termtype",
    "category": "Function",
    "text": "termtype(p::AbstractPolynomialLike)\n\nReturns the type of the terms of p.\n\ntermtype(::Type{PT}) where PT<:AbstractPolynomialLike\n\nReturns the type of the terms of a polynomial of type PT.\n\ntermtype(p::AbstractPolynomialLike, ::Type{T}) where T\n\nReturns the type of the terms of p but with coefficient type T.\n\ntermtype(::Type{PT}, ::Type{T}) where {PT<:AbstractPolynomialLike, T}\n\nReturns the type of the terms of a polynomial of type PT but with coefficient type T.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.coefficient",
    "page": "Reference",
    "title": "MultivariatePolynomials.coefficient",
    "category": "Function",
    "text": "coefficient(t::AbstractTermLike)\n\nReturns the coefficient of the term t.\n\nExamples\n\nCalling coefficient on 4x^2y should return 4.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.coefficienttype",
    "page": "Reference",
    "title": "MultivariatePolynomials.coefficienttype",
    "category": "Function",
    "text": "coefficient(p::AbstractPolynomialLike)\n\nReturns the coefficient type of p.\n\ncoefficient(::Type{PT}) where PT\n\nReturns the coefficient type of a polynomial of type PT.\n\nExamples\n\nCalling coefficienttype on (45)x^2y should return Rational{Int}, calling coefficienttype on 10x^2y + 20x should return Float64 and calling coefficienttype on xy should return Int.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.monomial",
    "page": "Reference",
    "title": "MultivariatePolynomials.monomial",
    "category": "Function",
    "text": "monomial(t::AbstractTermLike)\n\nReturns the monomial of the term t.\n\nExamples\n\nCalling monomial on 4x^2y should return x^2y.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.constantterm",
    "page": "Reference",
    "title": "MultivariatePolynomials.constantterm",
    "category": "Function",
    "text": "constantterm(α, p::AbstractPolynomialLike)\n\nCreates a constant term with coefficient α and the same variables as p.\n\nconstantterm(α, ::Type{PT} where {PT<:AbstractPolynomialType}\n\nCreates a constant term of the term type of a polynomial of type PT.\n\n\n\n"
},

{
    "location": "apireference.html#MultivariatePolynomials.zeroterm",
    "page": "Reference",
    "title": "MultivariatePolynomials.zeroterm",
    "category": "Function",
    "text": "zeroterm(p::AbstractPolynomialLike{T}) where T\n\nEquivalent to constantterm(zero(T), p).\n\nzeroterm(α, ::Type{PT} where {T, PT<:AbstractPolynomialType{T}}\n\nEquivalent to constantterm(zero(T), PT).\n\n\n\n"
},

{
    "location": "apireference.html#Terms-1",
    "page": "Reference",
    "title": "Terms",
    "category": "section",
    "text": "term\ntermtype\ncoefficient\ncoefficienttype\nmonomial\nconstantterm\nzeroterm"
},

{
    "location": "apireference.html#Polynomials-1",
    "page": "Reference",
    "title": "Polynomials",
    "category": "section",
    "text": "polynomial\npolynomialtype\nterms\nnterms\ncoefficients\nmonomials\nmindegree\nmaxdegree\nextdegree\nleadingterm\nleadingcoefficient\nleadingmonomial\nremoveleadingterm\nremovemonomials\nmonic"
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
    "location": "apireference.html#MultivariatePolynomials.emptymonovec",
    "page": "Reference",
    "title": "MultivariatePolynomials.emptymonovec",
    "category": "Function",
    "text": "emptymonovec(p::AbstractPolynomialLike)\n\nReturns an empty collection of the type of monomials(p).\n\n\n\n"
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
    "text": "monovec\nmonovectype\nemptymonovec\nsortmonovec\nmergemonovec"
},

]}

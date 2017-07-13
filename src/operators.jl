# We reverse the order of comparisons here so that the result
# of x < y is equal to the result of Monomial(x) < Monomial(y)
@pure isless(v1::AbstractVariable, v2::AbstractVariable) = name(v1) > name(v2)
isless(m1::AbstractTermLike, m2::AbstractTermLike) = isless(promote(m1, m2)...)

function isless(t1::AbstractTerm, t2::AbstractTerm)
    if monomial(t1) < monomial(t2)
        true
    elseif monomial(t1) == monomial(t2)
        coefficient(t1) < coefficient(t2)
    else
        false
    end
end

for op in [:+, :-, :*, :(==)]
    @eval $op(p1::APL, p2::APL) = $op(promote(p1, p2)...)
end
isapprox(p1::APL, p2::APL; kwargs...) = isapprox(promote(p1, p2)...; kwargs...)

# @eval $op(p::APL, α) = $op(promote(p, α)...) would be less efficient
for (op, fun) in [(:+, :plusconstant), (:-, :minusconstant), (:*, :multconstant), (:(==), :eqconstant)]
    @eval $op(p::APL, α) = $fun(p, α)
    @eval $op(α, p::APL) = $fun(α, p)
end
isapprox(p::APL, α; kwargs...) = isapproxconstant(promote(p, α)...; kwargs...)
isapprox(α, p::APL; kwargs...) = isapproxconstant(promote(p, α)...; kwargs...)

(-)(m::AbstractMonomialLike) = (-1) * m
(-)(t::AbstractTermLike) = (-coefficient(t)) * monomial(t)

# Avoid adding a zero constant that might artificially increase the Newton polytope
# Need to add polynomial conversion for type stability
plusconstant(p::APL, α) = iszero(α) ? polynomial(p) : p + term(α, p)
plusconstant(α, p::APL) = plusconstant(p, α)
minusconstant(p::APL, α) = iszero(α) ? polynomial(p) : p - term(α, p)
minusconstant(α, p::APL) = iszero(α) ? polynomial(-p) : term(α, p) - p

(+)(x::APL, y::MatPolynomial) = x + polynomial(y)
(+)(x::MatPolynomial, y::APL) = polynomial(x) + y
(+)(x::MatPolynomial, y::MatPolynomial) = polynomial(x) + polynomial(y)
(-)(x::APL, y::MatPolynomial) = x - polynomial(y)
(-)(x::MatPolynomial, y::APL) = polynomial(x) - y
(-)(x::MatPolynomial, y::MatPolynomial) = polynomial(x) - polynomial(y)

# Coefficients and variables commute
multconstant(α, v::AbstractVariable) = multconstant(α, monomial(v)) # TODO linear term
multconstant(m::AbstractMonomialLike, α) = multconstant(α, m)
multconstant(α, p::AbstractPolynomialLike) = multconstant(α, polynomial(p))
multconstant(p::AbstractPolynomialLike, α) = multconstant(polynomial(p), α)

multconstant(α, t::AbstractTerm)    = (α*coefficient(t)) * monomial(t)
multconstant(t::AbstractTerm, α)    = (coefficient(t)*α) * monomial(t)
(*)(m::AbstractMonomialLike, t::AbstractTerm) = coefficient(t) * (m * monomial(t))
(*)(t::AbstractTerm, m::AbstractMonomialLike) = coefficient(t) * (monomial(t) * m)
(*)(t1::AbstractTerm, t2::AbstractTerm) = (coefficient(t1) * coefficient(t2)) * (monomial(t1) * monomial(t2))

Base.transpose(v::AbstractVariable) = v
Base.transpose(m::AbstractMonomial) = m
Base.transpose(t::T) where {T <: AbstractTerm} = transpose(coefficient(t)) * monomial(t)
Base.transpose(p::AbstractPolynomialLike) = polynomial(map(transpose, terms(p)))

dot(p1::AbstractPolynomialLike, p2::AbstractPolynomialLike) = p1' * p2
dot(x, p::AbstractPolynomialLike) = x' * p
dot(p::AbstractPolynomialLike, x) = p' * x

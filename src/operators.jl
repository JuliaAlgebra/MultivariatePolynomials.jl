# We reverse the order of comparisons here so that the result
# of x < y is equal to the result of Monomial(x) < Monomial(y)
@pure isless(v1::AbstractVariable, v2::AbstractVariable) = name(v1) > name(v2)
isless(m1::AbstractTermLike, m2::AbstractTermLike) = isless(promote(m1, m2)...)

# Graded Lexicographic order
# First compare total degree, then lexicographic order
function isless(m1::AbstractMonomial{V}, m2::AbstractMonomial{V}) where {V}
    d1 = degree(m1)
    d2 = degree(m2)
    if d1 < d2
        return true
    elseif d1 > d2
        return false
    else
        return exponents(m1) < exponents(m2)
    end
end

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
    @eval $op(p1::AbstractPolynomialLike, p2::AbstractPolynomialLike) = $op(promote(p1, p2)...)
    @eval $op(p::AbstractPolynomialLike, x) = $op(promote(p, x)...)
    @eval $op(x, p::AbstractPolynomialLike) = $op(promote(x, p)...)
end

for op in [:+, :-]
    @eval $op(p1::AbstractTermLike, p2::AbstractTermLike) = $op(convert(AbstractPolynomial, p1), convert(AbstractPolynomial, p2))
    @eval $op(p1::AbstractPolynomial, p2::AbstractTermLike) = $op(p1, convert(AbstractPolynomial, p2))
    @eval $op(p1::AbstractTermLike, p2::AbstractPolynomial) = $op(convert(AbstractPolynomial, p1), p2)
end

(-)(t::AbstractTermLike) = -1 * t

(*)(t1::AbstractTerm, t2::AbstractTerm) = convert(AbstractTerm, coefficient(t1) * coefficient(t2), monomial(t1) * monomial(t2))
# TODO: this is inefficient
(*)(p1::AbstractPolynomial, p2::AbstractPolynomial) = sum(terms(p1) .* terms(p2).')
(*)(t::AbstractPolynomialLike) = t

@pure (==)(::AbstractVariable{N}, ::AbstractVariable{N}) where {N} = true
@pure (==)(::AbstractVariable, ::AbstractVariable) = false
(==)(m1::AbstractMonomial{V}, m2::AbstractMonomial{V}) where {V} = exponents(m1) == exponents(m2)
function (==)(t1::AbstractTerm, t2::AbstractTerm)
    c1 = coefficient(t1)
    c2 = coefficient(t2)
    if iszero(c1) && iszero(c2)
        true
    else
        c1 == c2 && monomial(t1) == monomial(t2)
    end
end
(==)(::AbstractPolynomialLike, ::Void) = false
(==)(::Void, ::AbstractPolynomialLike) = false

function compare_terms(p1::AbstractPolynomial, p2::AbstractPolynomial, op)
    i1 = 1
    i2 = 1
    t1 = terms(p1)
    t2 = terms(p2)
    while true
        while i1 <= length(t1) && coefficient(t1[i1]) == 0
            i1 += 1
        end
        while i2 <= length(t2) && coefficient(t2[i2]) == 0
            i2 += 1
        end
        if i1 > length(t1) && i2 > length(t2)
            return true
        end
        if i1 > length(t1) || i2 > length(t2)
            return false
        end
        if !op(t1[i1], t2[i2])
            return false
        end
        i1 += 1
        i2 += 1
    end
end

(==)(p1::AbstractPolynomial, p2::AbstractPolynomial) = compare_terms(p1, p2, (==))

isapprox(t1::AbstractTerm, t2::AbstractTerm; kwargs...) = isapprox(coefficient(t1), coefficient(t2); kwargs...) && monomial(t1) == monomial(t2)
isapprox(p1::AbstractPolynomial, p2::AbstractPolynomial; kwargs...) = compare_terms(p1, p2, (x, y) -> isapprox(x, y; kwargs...))

transpose(v::AbstractVariable) = v
transpose(m::AbstractMonomial) = m
transpose(t::T) where {T <: AbstractTerm} = T(coefficient(t)', monomial(t))
transpose(p::AbstractPolynomial) = convert(AbstractPolynomial, [transpose(t) for t in terms(p)])

dot(m::AbstractMonomialLike, x) = m * x
dot(x, m::AbstractMonomialLike) = x * m
dot(m1::AbstractMonomialLike, m2::AbstractMonomialLike) = m1 * m2

dot(t::AbstractTermLike, x) = dot(coefficient(t), x) * monomial(t)
dot(x, t::AbstractTermLike) = dot(x, coefficient(t)) * monomial(t)
dot(t1::AbstractTermLike, t2::AbstractTermLike) = dot(coefficient(t1), coefficient(t2)) * (monomial(t1) * monomial(t2))

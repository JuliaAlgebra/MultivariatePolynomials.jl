export variable

convertconstant(::Type{P}, α) where P<:APL = P(α)
Base.convert(::Type{P}, α) where P<:APL = convertconstant(P, α)
Base.convert(::Type{P}, p::APL) where P<:AbstractPolynomial = P(polynomial(p))

Base.convert(::Type{Any}, p::APL) = p
# Conversion polynomial -> scalar
function Base.convert(S::Type{<:Union{Number, T}}, p::APL{T}) where T
    s = zero(S)
    for t in terms(p)
        if !isconstant(t)
            # The polynomial is not constant
            throw(InexactError(:convert, S, p))
        end
        s += S(coefficient(t))
    end
    s
end

# Fix ambiguity caused by above conversions
Base.convert(::Type{P}, p::APL) where P<:APL = P(p)

Base.convert(::Type{PT}, p::PT) where {PT<:APL} = p
function Base.convert(::Type{MT}, t::AbstractTerm) where {MT<:AbstractMonomial}
    if isone(coefficient(t))
        monomial(t)
    else
        error("Cannot convert a term with a coefficient that is not one into a monomial")
    end
end

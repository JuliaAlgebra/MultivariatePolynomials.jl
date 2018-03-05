export variable

function convertconstant end
convert(::Type{P}, α) where P<:APL = convertconstant(P, α)
convert(::Type{P}, p::APL) where P<:AbstractPolynomial = convert(P, polynomial(p))

Base.convert(::Type{Any}, p::APL) = p
# Conversion polynomial -> scalar
function Base.convert(::Type{S}, p::APL) where {S}
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
Base.convert(::Type{PT}, p::PT) where {PT<:APL} = p
function Base.convert(::Type{MT}, t::AbstractTerm) where {MT<:AbstractMonomial}
    if isone(coefficient(t))
        monomial(t)
    else
        error("Cannot convert a term with a coefficient that is not one into a monomial")
    end
end

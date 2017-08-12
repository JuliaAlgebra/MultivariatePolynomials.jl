export variable

function convertconstant end
convert(::Type{P}, α) where P<:APL = convertconstant(P, α)
convert(::Type{P}, p::APL) where P<:AbstractPolynomial = convert(P, polynomial(p))

# Monomial -> variable
_errormono2var() = error("Monomial cannot be converted to a variable")
_mono2var() = _errormono2var()
function _checknovar() end
function _checknovar(ve, ves...)
    if iszero(ve[2])
        _checknovar(ves...)
    else
        _errormono2var()
    end
end
function _mono2var(ve, ves...)
    if iszero(ve[2])
        _mono2var(ves...)
    elseif isone(ve[2])
        _checknovar(ves...)
        ve[1]
    else
        _errormono2var()
    end
end
variable(m::AbstractMonomial) = _mono2var(powers(m)...)

Base.convert(::Type{Any}, p::APL) = p
# Conversion polynomial -> scalar
function Base.convert(::Type{S}, p::APL) where {S}
    s = zero(S)
    for t in terms(p)
        if !isconstant(t)
            # The polynomial is not constant
            throw(InexactError())
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


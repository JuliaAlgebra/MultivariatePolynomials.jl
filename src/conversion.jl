export variable

function convertconstant end
Base.convert(::Type{P}, α) where P<:APL = convertconstant(P, α)
function Base.convert(::Type{P}, p::P) where {T, P<:AbstractPolynomial{T}}
    return p
end
function Base.convert(::Type{P}, p::APL) where {T, P<:AbstractPolynomial{T}}
    return convert(P, polynomial(p, T))
end

function Base.convert(::Type{V}, mono::AbstractMonomial) where V <: AbstractVariable
    variable = nothing
    for v in variables(mono)
        d = degree(mono, v)
        if isone(d)
            if variable === nothing
                variable = v
            else
                throw(InexactError(:convert, V, mono))
            end
        elseif !iszero(d)
            throw(InexactError(:convert, V, mono))
        end
    end
    if variable === nothing
        throw(InexactError(:convert, V, mono))
    end
    return variable
end

function Base.convert(::Type{M}, t::AbstractTerm) where M <: AbstractMonomialLike
    if isone(coefficient(t))
        return convert(M, monomial(t))
    else
        throw(InexactError(:convert, M, t))
    end
end
function Base.convert(TT::Type{<:AbstractTerm{T}}, m::AbstractMonomialLike) where T
    return convert(TT, one(T) * m)
end
function Base.convert(TT::Type{<:AbstractTerm{T}}, t::AbstractTerm) where T
    return convert(TT, convert(T, coefficient(t)) * monomial(t))
end

# Base.convert(::Type{T}, t::T) where {T <: AbstractTerm} is ambiguous with above method.
# we need the following:
function Base.convert(::Type{TT}, t::TT) where {T, TT <: AbstractTerm{T}}
    return t
end

function Base.convert(::Type{T}, p::AbstractPolynomial) where T <: AbstractTermLike
    if iszero(nterms(p))
        convert(T, zeroterm(p))
    elseif isone(nterms(p))
        convert(T, leadingterm(p))
    else
        throw(InexactError(:convert, T, p))
    end
end

MA.scaling(p::AbstractPolynomialLike{T}) where {T} = convert(T, p)
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

Base.convert(::Type{PT}, p::PT) where {PT<:APL} = p

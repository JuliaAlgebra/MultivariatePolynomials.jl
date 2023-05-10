export variable, convert_to_constant

function convert_constant end
Base.convert(::Type{P}, α) where {P<:APL} = convert_constant(P, α)
function convert_constant(::Type{TT}, α) where {T,TT<:AbstractTerm{T}}
    return term(convert(T, α), constant_monomial(TT))
end
function convert_constant(::Type{PT}, α) where {PT<:AbstractPolynomial}
    return convert(PT, convert(term_type(PT), α))
end
function Base.convert(::Type{P}, p::APL) where {T,P<:AbstractPolynomial{T}}
    return error("`convert` not implemented for $P")
end

function Base.convert(
    ::Type{V},
    mono::AbstractMonomial,
) where {V<:AbstractVariable}
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

function Base.convert(
    ::Type{M},
    t::AbstractTerm,
) where {M<:AbstractMonomialLike}
    if isone(coefficient(t))
        return convert(M, monomial(t))
    else
        throw(InexactError(:convert, M, t))
    end
end
function Base.convert(
    TT::Type{<:AbstractTerm{T}},
    m::AbstractMonomialLike,
) where {T}
    return convert(TT, term(one(T), convert(monomial_type(TT), m)))
end
function Base.convert(TT::Type{<:AbstractTerm{T}}, t::AbstractTerm) where {T}
    return convert(
        TT,
        term(
            convert(T, coefficient(t)),
            convert(monomial_type(TT), monomial(t)),
        ),
    )
end

# Base.convert(::Type{T}, t::T) where {T <: AbstractTerm} is ambiguous with above method.
# we need the following:
function Base.convert(::Type{TT}, t::TT) where {T,TT<:AbstractTerm{T}}
    return t
end

function Base.convert(
    ::Type{T},
    p::AbstractPolynomial,
) where {T<:AbstractTermLike}
    if iszero(nterms(p))
        convert(T, zero_term(p))
    elseif isone(nterms(p))
        convert(T, leading_term(p))
    else
        throw(InexactError(:convert, T, p))
    end
end

MA.scaling(p::AbstractPolynomialLike{T}) where {T} = convert(T, p)
# Conversion polynomial -> constant
# We don't define a method for `Base.convert` to reduce invalidations;
# see https://github.com/JuliaAlgebra/MultivariatePolynomials.jl/pull/172
function convert_to_constant(::Type{S}, p::APL) where {S}
    s = zero(S)
    for t in terms(p)
        if !isconstant(t)
            # The polynomial is not constant
            throw(InexactError(:convert_to_constant, S, p))
        end
        s = MA.add!!(s, convert(S, coefficient(t)))
    end
    return s
end
Base.convert(::Type{T}, p::APL) where {T<:Number} = convert_to_constant(T, p)
function convert_to_constant(p::APL{S}) where {S}
    return convert_to_constant(S, p)
end

# Also covers, e.g., `convert(APL, ::P)` where `P<:APL`
Base.convert(::Type{PT}, p::PT) where {PT<:APL} = p

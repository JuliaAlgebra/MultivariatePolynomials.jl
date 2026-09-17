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
    t::SA.Term,
) where {M<:AbstractMonomialLike}
    if isone(coefficient(t))
        return convert(M, monomial(t))
    else
        throw(InexactError(:convert, M, t))
    end
end

function Base.convert(
    ::Type{T},
    p::AbstractPolynomial,
) where {T<:AbstractMonomialLike}
    if iszero(nterms(p))
        convert(T, zero_term(p))
    elseif isone(nterms(p))
        convert(T, leading_term(p))
    else
        throw(InexactError(:convert, T, p))
    end
end
MA.scaling(p::AbstractPolynomialLike) = convert(coefficient_type(p), p)
# Conversion polynomial -> constant
# We don't define a method for `Base.convert` to reduce invalidations;
# see https://github.com/JuliaAlgebra/MultivariatePolynomials.jl/pull/172
function convert_to_constant(::Type{S}, p::_APL) where {S}
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
Base.convert(::Type{T}, p::_APL) where {T<:Number} = convert_to_constant(T, p)
function convert_to_constant(p::_APL)
    return convert_to_constant(coefficient_type(p), p)
end

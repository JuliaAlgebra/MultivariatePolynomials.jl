# Base on the TypedPolynomials/abstract/substition.jl written by Robin Deits

# TODO Vararg{<:...} -> Vararg{...}
const Substitution = Pair{<:AbstractVariable}
const MultiSubstitution{N} =
    Pair{<:Tuple{Vararg{AbstractVariable,N}},<:Tuple{Vararg{Any,N}}}
const MultiVectorSubstitution =
    Pair{<:Tuple{Vararg{AbstractVariable}},<:AbstractVector}

# When the variables are promoted to be in the same vector they could be promoted into a monomial
const VectorMultiSubstitution =
    Pair{<:AbstractVector{<:AbstractMonomialLike},<:Tuple}
const VectorMultiVectorSubstitution =
    Pair{<:AbstractVector{<:AbstractMonomialLike},<:AbstractVector}
const AbstractMultiSubstitution = Union{
    MultiSubstitution,
    MultiVectorSubstitution,
    VectorMultiVectorSubstitution,
    VectorMultiSubstitution,
}
const AbstractSubstitution = Union{Substitution,AbstractMultiSubstitution}
const Substitutions = Tuple{Vararg{AbstractSubstitution}}

abstract type AbstractSubstitutionType end
struct Subs <: AbstractSubstitutionType end
struct Eval <: AbstractSubstitutionType end
const _AST = AbstractSubstitutionType

"""
    subs(polynomial, (x, y)=>(1, 2))

is equivalent to:

    subs(polynomial, (x=>1, y=>2))
"""
function substitute(st::_AST, p::_APL, s::AbstractMultiSubstitution)
    return substitute(st, p, (s,))
end

_replacement_type(::Type{Pair{V,T}}) where {V<:AbstractVariable,T} = T
function _replacement_type(::Type{Pair{V,T}}) where {V,T}
    return _replacement_eltype(T)
end
_replacement_eltype(::Type{T}) where {T} = eltype(T)
function _replacement_eltype(::Type{T}) where {T<:Tuple}
    return promote_type(fieldtypes(T)...)
end

function _substitution_type(
    ::Type{Eval},
    ::Type,
    args::Vararg{Type,N},
) where {N}
    return N == 0 ? Int : promote_type(map(_replacement_type, args)...)
end
function _substitution_type(
    ::Type{Subs},
    ::Type{P},
    args::Vararg{Type,N},
) where {P,N}
    return promote_type(variable_union_type(P), map(_replacement_type, args)...)
end

function _substitution_value(::Type{T}, ::Subs, v) where {T}
    return convert(T, v)
end
function _substitution_value(::Type, ::Eval, v)
    return throw(ArgumentError("No value provided for variable $v"))
end

function _substitution_matches(v, w)
    return v == w ||
           (w isa AbstractVariable && !isreal(w) && ordinary_variable(v) == w)
end

function _substitution_part(v, w, value)
    v == w && return value
    # A value for a complex variable also determines its conjugate and parts.
    isconj(v) && return conj(value)
    isrealpart(v) && return real(value)
    return imag(value)
end

function _substitution_value(
    ::Type{T},
    st::_AST,
    v,
    s::Substitution,
    tail::Vararg{AbstractSubstitution,N},
) where {T,N}
    if _substitution_matches(v, first(s))
        return convert(T, _substitution_part(v, first(s), last(s)))
    end
    return _substitution_value(T, st, v, tail...)
end
function _substitution_value(
    ::Type{T},
    st::_AST,
    v,
    s::AbstractMultiSubstitution,
    tail::Vararg{AbstractSubstitution,N},
) where {T,N}
    i = findfirst(Base.Fix1(_substitution_matches, v), first(s))
    if !isnothing(i)
        return convert(T, _substitution_part(v, first(s)[i], last(s)[i]))
    end
    return _substitution_value(T, st, v, tail...)
end

_check_substitution(::Substitution) = nothing
function _check_substitution(s::AbstractMultiSubstitution)
    length(first(s)) == length(last(s)) || throw(
        ArgumentError(
            "Substitution variables and values must have equal lengths",
        ),
    )
    return
end

function substitute(st::_AST, p::_APL, s::Substitutions)
    foreach(_check_substitution, s)
    if st isa Eval && isempty(s)
        return convert_to_constant(p)
    end
    T =
        variables(p) isa Tuple && st isa Subs ? Any :
        _substitution_type(typeof(st), typeof(p), typeof.(s)...)
    if variables(p) isa Tuple
        vals = map(v -> _substitution_value(T, st, v, s...), variables(p))
        return _substitute(p, vals)
    end
    # Prepare homogeneous values once, also when p has no nonzero terms.
    vals = T[_substitution_value(T, st, v, s...) for v in variables(p)]
    return _substitute(p, vals)
end
function substitute(st::_AST, p::_APL, s::AbstractSubstitution...)
    return substitute(st, p, s)
end

_substitute(::AbstractVariable, vals) = only(vals)

_substitute_powers(::Tuple{}, ::Tuple{}) = 1
function _substitute_powers(vals::Tuple{Any}, exps::Tuple{Int})
    return only(vals)^only(exps)
end
function _substitute_powers(vals::Tuple, exps::Tuple)
    return first(vals)^first(exps) *
           _substitute_powers(Base.tail(vals), Base.tail(exps))
end
function _substitute(m::AbstractMonomial, vals::Tuple)
    return _substitute_powers(vals, exponents(m))
end

function _substitute(m::AbstractMonomial, vals::AbstractVector{T}) where {T}
    if isempty(vals)
        return one(MA.promote_operation(*, T, T))
    end
    z = exponents(m)
    result = vals[1]^z[1]
    for i in 2:length(vals)
        # Keep zero powers: their parents are part of the result's basis.
        result = result * vals[i]^z[i]
    end
    return result
end

function _substitute(t::AbstractTerm, vals)
    return coefficient(t) * _substitute(monomial(t), vals)
end

function MA.promote_operation(
    ::typeof(substitute),
    ::Type{S},
    ::Type{M},
    args::Vararg{Type,N},
) where {S<:_AST,M<:AbstractMonomialLike,N}
    T = _substitution_type(S, M, args...)
    return M <: AbstractVariable ? T : MA.promote_operation(*, T, T)
end

function MA.promote_operation(
    ::typeof(substitute),
    ::Type{S},
    ::Type{T},
    args::Vararg{Type,N},
) where {S<:AbstractSubstitutionType,T<:AbstractTerm,N}
    M = MA.promote_operation(substitute, S, monomial_type(T), args...)
    U = coefficient_type(T)
    return MA.promote_operation(*, U, M)
end

## Polynomials
_polynomial(α) = α
_polynomial(p::_APL) = polynomial(p)
function _substitute(p::AbstractPolynomial, vals)
    # Evaluating the constant monomial prepares the result's common basis,
    # including for zero polynomials. Subsequent additions can mutate it.
    # Noncommutative monomial products may drop zero-degree variables. Algebra
    # elements keep their parents when raised to the zeroth power.
    unit_vals = is_commutative(p) ? vals : _polynomial.(vals)
    unit = _substitute(constant_monomial(p), unit_vals)
    result = _polynomial(zero(coefficient_type(p)) * unit)
    for t in terms(p)
        # Noncommutative keys have their own variable sequence, with repetitions.
        term_vals =
            is_commutative(t) ? vals :
            [vals[findfirst(isequal(v), variables(p))] for v in variables(t)]
        result = MA.operate!!(+, result, _substitute(t, term_vals))
    end
    return result
end

function MA.promote_operation(
    ::typeof(substitute),
    ::Type{S},
    ::Type{P},
    args::Vararg{Type,N},
) where {S<:AbstractSubstitutionType,P<:AbstractPolynomial,N}
    T = MA.promote_operation(substitute, S, term_type(P), args...)
    return MA.promote_operation(+, T, T)
end

## Fallbacks
function substitute(st::_AST, q::RationalPoly, s::Substitutions)
    return substitute(st, q.num, s) / substitute(st, q.den, s)
end

# Substitute Arrays
function substitute(st::_AST, A::AbstractArray{<:_APL}, s::Substitutions)
    return map(p -> substitute(st, p, s), A)
end
## Everything else
substitute(::_AST, x, s::Substitutions) = x

"""
    subs(p, s::AbstractSubstitution...)

Apply the substitutions `s` to `p`.
Use `p(s...)` if we are sure that all the variables are substited in `s`.

The allowed substutions are:

* `v => p` where `v` is a variable and `p` a polynomial, e.g. `x => 1` or `x => x^2*y + x + y`.
* `V => P` where `V` is a tuple or vector of variables and `P` a tuple or vector of polynomials, e.g. `(x, y) => (y, x)` or `(y, x) => (2, 1)`.

The order of the variables is lexicographic with the name with TypedPolynomials and by order of creation with DynamicPolynomials.
Since there is no guarantee on the order of the variables, substitution directly with a tuple or a vector is not allowed.
You can use `p(variables(p) => (1, 2))` instead if you are sure of the order of the variables (e.g. the name order matches the creation order).

### Examples

```julia
p = 3x^2*y + x + 2y + 1
p(x => 2, y => 1) # Return type is Int
subs(p, x => 2, y => 1) # Return type is Int in TypedPolynomials but is a polynomial of Int coefficients in DynamicPolynomials
subs(p, y => x*y^2 + 1)
p(y => 2) # Use subs(p, y => 2) for partial substitution.
```

"""
subs(p, s::AbstractSubstitution...) = substitute(Subs(), p, s)

(p::RationalPoly)(s::AbstractSubstitution...) = p.num(s...) / p.den(s...)

function (p::AbstractPolynomialLike)(s::AbstractSubstitution...)
    return substitute(Eval(), p, s)
end

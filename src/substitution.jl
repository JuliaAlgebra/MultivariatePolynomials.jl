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
function _monomial_vector_to_variable_tuple(
    s::Pair{<:AbstractVector{<:AbstractMonomial}},
)
    return variable.(Tuple(s.first)) => s.second
end
_monomial_vector_to_variable_tuple(s) = s

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
    return substitute(st, p, _flatten_subs(s))
end

## Variables
function substitute(st::_AST, v::AbstractVariable, s::Substitutions)
    return substitute(st, v, s...)
end

## Monomials
function powersubstitute(
    st::_AST,
    s::Substitutions,
    p::Tuple{AbstractVariable,Integer},
)
    return substitute(st, p[1], s...)^p[2]
end
function powersubstitute(
    st::_AST,
    s::Substitutions,
    p::Tuple{AbstractVariable,Integer},
    p2...,
)
    return powersubstitute(st, s, p) * powersubstitute(st, s, p2...)
end

function _promote_subs(
    S,
    ::Type{V},
    s::Substitution,
) where {V<:AbstractVariable}
    return MA.promote_operation(substitute, S, V, typeof(s))
end

function _flatten_subs(s::AbstractMultiSubstitution)
    return pair_zip(_monomial_vector_to_variable_tuple(s))
end

function _flatten_subs(s::Substitution)
    return (s,)
end

# Turn a tuple of `AbstractSubstitution` into a `Tuple` if `Substitution`
function _flatten_subs(
    s::AbstractSubstitution,
    tail::Vararg{AbstractSubstitution,N},
) where {N}
    return (_flatten_subs(s)..., _flatten_subs(tail...)...)
end

function power_promote(
    S,
    ::Type{V},
    s::Substitutions,
) where {V<:AbstractVariable}
    T = MA.promote_operation(substitute, S, V, typeof.(_flatten_subs(s...))...)
    return MA.promote_operation(*, T, T)
end

function power_promote(
    S,
    ::Vector{V},
    s::Substitutions,
) where {V<:AbstractVariable}
    return _promote_subs_power(S, V, s)
end

function power_promote(
    S,
    ::Tuple{V},
    s::Substitutions,
) where {V<:AbstractVariable}
    return _promote_subs_power(S, V, s)
end

function power_promote(
    S,
    vars::Tuple{V,Vararg{AbstractVariable,N}},
    s::Substitutions,
) where {V<:AbstractVariable,N}
    return MA.promote_operation(
        *,
        _promote_subs_power(S, V, s),
        _promote_subs_mono(S, Base.tail(vars), s),
    )
end

function substitute(st::_AST, m::AbstractMonomial, s::Substitutions)
    if isconstant(m)
        return one(power_promote(typeof(st), variables(m), s))
    else
        return powersubstitute(st, s, powers(m)...)
    end
end

## Terms
function substitute(st::_AST, t::AbstractTerm, s::Substitutions)
    return coefficient(t) * substitute(st, monomial(t), s)
end

function MA.promote_operation(
    ::typeof(substitute),
    ::Type{Eval},
    ::Type{M},
    ::Type{Pair{V,T}},
) where {M<:AbstractMonomial,V<:AbstractVariable,T}
    return MA.promote_operation(*, T, T)
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
function substitute(st::_AST, p::AbstractPolynomial, s::Substitutions)
    if iszero(p)
        _polynomial(substitute(st, zero_term(p), s))
    else
        ts = terms(p)
        r1 = substitute(st, ts[1], s)
        R = MA.promote_operation(+, typeof(r1), typeof(r1))
        result::R = convert(R, r1)
        for i in 2:length(ts)
            result += substitute(st, ts[i], s)
        end
        result
    end
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
function substitute(st::_AST, p::_APL, s::Substitutions)
    return substitute(st, polynomial(p), s)
end
function substitute(st::_AST, q::RationalPoly, s::Substitutions)
    return substitute(st, q.num, s) / substitute(st, q.den, s)
end

# subs(x, x=>x+y, y=>2) would call substitute(Subs(), x+y, y=>2)
#substitute(st::_AST, p::Union{_APL, RationalPoly}, s::AbstractSubstitution...) = substitute(st, p, s)

# Substitute Arrays
function substitute(st::_AST, A::AbstractArray{<:_APL}, s::Substitutions)
    return map(p -> substitute(st, p, s), A)
end
## Everything else
substitute(::_AST, x, s::Substitutions) = x
# subs(x, x=>1, y=>2) would call substitute(Subs(), 1, y=>2)
#substitute(::_AST, x, s::AbstractSubstitution...) = x

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
p(y => 2) # Do not do that, this works fine with TypedPolynomials but it will not return a correct result with DynamicPolynomials since it thinks that the return type is `Int`.
```

"""
subs(p, s::AbstractSubstitution...) = substitute(Subs(), p, s)

(p::RationalPoly)(s::AbstractSubstitution...) = p.num(s...) / p.den(s...)

# Minimal implementation of commutative variables so that the tests of
# MultivariatePolynomials do not depend on DynamicPolynomials or
# TypedPolynomials. Terms and polynomials are the ones of StarAlgebras.
module SimplePolynomials

import MultivariatePolynomials as MP
import MutableArithmetics as MA
import StarAlgebras as SA

export Variable, @polyvar, @complex_polyvar

@enum ComplexKind REAL COMPLEX CONJ REAL_PART IMAG_PART

const _COUNTER = Ref(0)

# Variables created first are larger, the variables of a monomial are sorted
# in decreasing order.
struct Variable <: MP.AbstractVariable
    name::String
    id::Int
    kind::ComplexKind
end
function Variable(name::AbstractString, kind::ComplexKind = REAL)
    return Variable(name, _COUNTER[] += 1, kind)
end

function _array_variable(prefix, kind, indices...)
    return map(Iterators.product(indices...)) do i
        return Variable("$prefix[$(join(i, ","))]", kind)
    end
end

function _assign(var::Symbol, kind)
    return var, :($(esc(var)) = $Variable($(string(var)), $kind))
end
function _assign(var::Expr, kind)
    Base.Meta.isexpr(var, :ref) ||
        error("Expected $var to be of the form varname[idxset]")
    name = var.args[1]
    indices = esc.(var.args[2:end])
    return name,
    :($(esc(name)) = $_array_variable($(string(name)), $kind, $(indices...)))
end

function _polyvar(args, kind)
    vars_exprs = _assign.(args, kind)
    return Expr(
        :block,
        last.(vars_exprs)...,
        Expr(:tuple, esc.(first.(vars_exprs))...),
    )
end

macro polyvar(args...)
    return _polyvar(args, REAL)
end

macro complex_polyvar(args...)
    return _polyvar(args, COMPLEX)
end

Base.hash(x::Variable, u::UInt) = hash(x.id, u)
Base.broadcastable(x::Variable) = Ref(x)

function Base.cmp(x::Variable, y::Variable)
    if x.id == y.id
        return cmp(y.kind, x.kind)
    end
    return cmp(y.id, x.id)
end
Base.:(==)(x::Variable, y::Variable) = iszero(cmp(x, y))

MP.name(v::Variable) = v.name
function MP.name_base_indices(v::Variable)
    splits = split(v.name, r"[\[,\]]\s*", keepempty = false)
    if length(splits) == 1
        return v.name, Int[]
    end
    return splits[1], parse.(Int, splits[2:end])
end

MP.variables(v::Variable) = [v]
MP.exponents(::Variable) = [1]
MP.ordering(::Union{Variable,Type{Variable}}) = MP.Graded{MP.LexOrder}()
MP.is_commutative(::Type{<:Union{Variable,Vector{Variable}}}) = true

Base.isreal(x::Variable) = x.kind != COMPLEX && x.kind != CONJ
MP.isrealpart(x::Variable) = x.kind == REAL_PART
MP.isimagpart(x::Variable) = x.kind == IMAG_PART
MP.isconj(x::Variable) = x.kind == CONJ
function MP.ordinary_variable(x::Variable)
    return x.kind == REAL || x.kind == COMPLEX ? x :
           Variable(x.name, x.id, COMPLEX)
end
function Base.conj(x::Variable)
    x.kind == COMPLEX && return Variable(x.name, x.id, CONJ)
    x.kind == CONJ && return Variable(x.name, x.id, COMPLEX)
    return x
end
function Base.real(x::Variable)
    isreal(x) && return x
    return Variable(x.name, x.id, REAL_PART)
end
function Base.imag(x::Variable)
    x.kind == COMPLEX && return Variable(x.name, x.id, IMAG_PART)
    x.kind == CONJ &&
        return SA.Term(-1, MP.monomial(Variable(x.name, x.id, IMAG_PART)))
    return MA.Zero()
end

const Monomial = MP.BasisElement{MP.Monomial,Vector{Variable},Vector{Int}}

function MP.variable_union_type(
    ::Union{M,Type{M}} where {M<:Union{Variable,Monomial}},
)
    return Variable
end
MP.monomial_type(::Union{Variable,Type{Variable}}) = Monomial
MP.monomial(v::Variable) = MP.BasisElement{MP.Monomial}(v)
function MP.monomial(vars::Vector{Variable}, z::Vector{Int})
    length(vars) == length(z) ||
        throw(ArgumentError("There should be as many exponents as variables"))
    @assert issorted(vars, rev = true)
    return MP.BasisElement(MP.Variables{MP.Monomial}(vars), z)
end
function MP.constant_monomial(::Type{Monomial})
    return MP.monomial(Variable[], Int[])
end
Base.convert(::Type{Monomial}, x::Variable) = MP.monomial(x)
Base.:^(x::Variable, i::Int) = MP.monomial([x], [i])
Base.:+(x::Variable, y::Variable) = MP.term(x) + MP.term(y)
Base.:-(x::Variable, y::Variable) = MP.term(x) - MP.term(y)

function MP.similar_variable(
    ::Union{Variable,Type{Variable}},
    ::Type{Val{S}},
) where {S}
    return Variable(string(S))
end

function SA.promote_with_map(
    ::Variable,
    all_vars::Vector{Variable},
    map::MP.ExponentMap,
)
    return MP.monomial(copy(all_vars), map([1])), map
end

function MP.divides(m1::Monomial, m2::Monomial)
    v2 = MP.variables(m2)
    e2 = MP.exponents(m2)
    for (v, e) in zip(MP.variables(m1), MP.exponents(m1))
        iszero(e) && continue
        i = findfirst(isequal(v), v2)
        (isnothing(i) || e2[i] < e) && return false
    end
    return true
end

const _Algebra = typeof(MP.algebra(MP.FullBasis{MP.Monomial}(Variable[])))
const Polynomial{T} = SA.AlgebraElement{
    T,
    _Algebra,
    SA.SparseCoefficients{
        Vector{Int},
        T,
        Vector{Vector{Int}},
        Vector{T},
        MP.Graded{MP.LexOrder},
    },
}
const Term{T} = SA.Term{T,_Algebra,Vector{Int}}

end

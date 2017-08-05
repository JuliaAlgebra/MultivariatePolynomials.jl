export subs

# TODO Vararg{<:...} -> Vararg{...}
const Substitution = Pair{<:AbstractVariable}
const MultiSubstitution{N} = Pair{<:Tuple{Vararg{AbstractVariable, N}}, <:Tuple{Vararg{<:Any, N}}}
const MultiVectorSubstitution = Pair{<:Tuple{Vararg{AbstractVariable}}, <:AbstractVector}
const VectorMultiSubstitution = Pair{<:AbstractVector{<:AbstractVariable}, <:Tuple}
const VectorMultiVectorSubstitution = Pair{<:AbstractVector{<:AbstractVariable}, <:AbstractVector}

const AbstractMultiSubstitution = Union{MultiSubstitution, MultiVectorSubstitution, VectorMultiVectorSubstitution, VectorMultiSubstitution}
const AbstractSubstitution = Union{Substitution, AbstractMultiSubstitution}
const Substitutions = Tuple{Vararg{AbstractSubstitution}}

abstract type AbstractSubstitutionType end
struct Subs <: AbstractSubstitutionType end
struct Eval <: AbstractSubstitutionType end
const AST = AbstractSubstitutionType

## Monomials
substitute(st::AST, s::Substitutions, p::Tuple{AbstractVariable, Integer}) = substitute(st, p[1], s...)^p[2]
substitute(st::AST, s::Substitutions, p::Tuple{AbstractVariable, Integer}, p2...) = substitute(st, s, p) * substitute(st, s, p2...)
substitute(st::AST, m::AbstractMonomial, s::Substitutions) = substitute(st, s, powers(m)...)

## Terms
substitute(st::AST, t::AbstractTerm, s::Substitutions) = coefficient(t) * substitute(st, monomial(t), s)

## Polynomials
function substitute(st::AST, p::AbstractPolynomial, s::Substitutions)
    ts = terms(p)
    @assert length(ts) >= 1
    r1 = substitute(st, ts[1], s)
    R = Base.promote_op(+, typeof(r1), typeof(r1))
    result::R = convert(R, r1)
    for i in 2:length(ts)
        result += substitute(st, ts[i], s)
    end
    result
end

## Varargs catchers
substitute(st::AST, p::APL, s::Substitutions) = substitute(st, polynomial(p), s)
substitute(st::AST, q::RationalPoly, s::Substitutions) = substitute(st, q.num, s) / substitute(st, q.den, s)

## Everything else
substitute(::AST, x, s::Substitutions) = x

"""
    subs(polynomial, (x, y)=>(1, 2))

is equivalent to:

    subs(polynomial, (x=>1, y=>2))
"""
substitute(::AST, p, s::MultiSubstitution) = subs(p, pairzip(s))

# inefficient but convenient method to allow subs(p, (x, y)=>[1, 2])
substitute(::AST, p, s::MultiVectorSubstitution) = subs(p, s.first => Tuple(s.second))

"""
eval replace all variables by a new value and subs replace some of them.
we use different function for that because of type inferability.
with eval, Julia knows that all PolyVar will be replaced by values so it can do better inference.
"""
subs(p, s::AbstractSubstitution...) = substitute(Subs(), p, s)

(p::MatPolynomial)(s::AbstractSubstitution...) = polynomial(p)(s...)
(p::RationalPoly)(s::AbstractSubstitution...) = p.num(s...) / p.den(s...)

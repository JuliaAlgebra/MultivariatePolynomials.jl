# Base on the TypedPolynomials/abstract/substition.jl written by Robin Deits
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

"""
    subs(polynomial, (x, y)=>(1, 2))

is equivalent to:

    subs(polynomial, (x=>1, y=>2))
"""
substitute(st::AST, p::APL, s::MultiSubstitution) = substitute(st, p, pairzip(s))

# inefficient but convenient method to allow subs(p, (x, y)=>[1, 2])
substitute(st::AST, p::APL, s::MultiVectorSubstitution) = substitute(st, p, s.first => Tuple(s.second))
substitute(st::AST, p::APL, s::VectorMultiSubstitution) = substitute(st, p, Tuple(s.first) => s.second)
substitute(st::AST, p::APL, s::VectorMultiVectorSubstitution) = substitute(st, p, Tuple(s.first) => Tuple(s.second))


# Evaluate the stream
# If I do s2..., then
# subs(x, x=>x+y, y=>2) would call substitute(Subs(), x+y, y=>2)
# subs(x, x=>1, y=>2) would call substitute(Subs(), 1, y=>2)
# so it would force use to also define
# sustitute(::AST, ..., ::AbstractSubstitution...) for Any, APL and RationalPoly.
#substitute(st::AST, p::APL, s1::AbstractSubstitution, s2::AbstractSubstitution...) = substitute(st, substitute(st, p, s1), s2...)
substitute(st::AST, p::APL, s1::AbstractSubstitution, s2::AbstractSubstitution...) = substitute(st, substitute(st, p, s1), s2)

## Variables
substitute(st::AST, v::AbstractVariable, s::Substitutions) = substitute(st, v, s...)

## Monomials
powersubstitute(st::AST, s::Substitutions, p::Tuple{AbstractVariable, Integer}) = substitute(st, p[1], s...)^p[2]
powersubstitute(st::AST, s::Substitutions, p::Tuple{AbstractVariable, Integer}, p2...) = powersubstitute(st, s, p) * powersubstitute(st, s, p2...)
substitute(st::AST, m::AbstractMonomial, s::Substitutions) = powersubstitute(st, s, powers(m)...)

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

## Fallbacks
substitute(st::AST, p::APL, s::Substitutions) = substitute(st, polynomial(p), s)
substitute(st::AST, q::RationalPoly, s::Substitutions) = substitute(st, q.num, s) / substitute(st, q.den, s)

# subs(x, x=>x+y, y=>2) would call substitute(Subs(), x+y, y=>2)
#substitute(st::AST, p::Union{APL, RationalPoly}, s::AbstractSubstitution...) = substitute(st, p, s)

## Everything else
substitute(::AST, x, s::Substitutions) = x
# subs(x, x=>1, y=>2) would call substitute(Subs(), 1, y=>2)
#substitute(::AST, x, s::AbstractSubstitution...) = x

"""
eval replace all variables by a new value and subs replace some of them.
we use different function for that because of type inferability.
with eval, Julia knows that all PolyVar will be replaced by values so it can do better inference.
"""
subs(p, s::AbstractSubstitution...) = substitute(Subs(), p, s)

(p::MatPolynomial)(s::AbstractSubstitution...) = polynomial(p)(s...)
(p::RationalPoly)(s::AbstractSubstitution...) = p.num(s...) / p.den(s...)

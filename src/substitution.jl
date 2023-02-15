# Base on the TypedPolynomials/abstract/substition.jl written by Robin Deits
export subs

# TODO Vararg{<:...} -> Vararg{...}
const Substitution = Pair{<:AbstractVariable}
const MultiSubstitution{N} = Pair{<:Tuple{Vararg{AbstractVariable, N}}, <:Tuple{Vararg{Any, N}}}
const MultiVectorSubstitution = Pair{<:Tuple{Vararg{AbstractVariable}}, <:AbstractVector}

# When the variables are promoted to be in the same vector they could be promoted into a monomial
const VectorMultiSubstitution = Pair{<:AbstractVector{<:AbstractMonomialLike}, <:Tuple}
const VectorMultiVectorSubstitution = Pair{<:AbstractVector{<:AbstractMonomialLike}, <:AbstractVector}
_monov2vart(s::Pair{<:AbstractVector{<:AbstractMonomial}}) = variable.(Tuple(s.first)) => s.second
_monov2vart(s) = s

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
substitute(st::AST, p::APL, s::AbstractMultiSubstitution) = substitute(st, p, pairzip(_monov2vart(s)))

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
_polynomial(α) = α
_polynomial(p::APL) = polynomial(p)
function substitute(st::AST, p::AbstractPolynomial, s::Substitutions)
    if iszero(p)
        _polynomial(substitute(st, zeroterm(p), s))
    else
        ts = terms(p)
        r1 = substitute(st, ts[1], s)
        R = Base.promote_op(+, typeof(r1), typeof(r1))
        result::R = convert(R, r1)
        for i in 2:length(ts)
            result += substitute(st, ts[i], s)
        end
        result
    end
end

## Fallbacks
substitute(st::AST, p::APL, s::Substitutions) = substitute(st, polynomial(p), s)
substitute(st::AST, q::RationalPoly, s::Substitutions) = substitute(st, q.num, s) / substitute(st, q.den, s)

# subs(x, x=>x+y, y=>2) would call substitute(Subs(), x+y, y=>2)
#substitute(st::AST, p::Union{APL, RationalPoly}, s::AbstractSubstitution...) = substitute(st, p, s)

# Substitute Arrays
function substitute(st::AST, A::AbstractArray{<:APL}, s::Substitutions)
    map(p -> substitute(st, p, s), A)
end
## Everything else
substitute(::AST, x, s::Substitutions) = x
# subs(x, x=>1, y=>2) would call substitute(Subs(), 1, y=>2)
#substitute(::AST, x, s::AbstractSubstitution...) = x

"""
    subs(p, s::AbstractSubstitution...)

Apply the substitutions `s` to `p`.
Use `p(s...)` if we are sure that all the variables are substited in `s`.

The allowed substutions are:

* `v => p` where `v` is a variable and `p` a polynomial, e.g. `x => 1` or `x => x^2*y + x + y`.
* `V => P` where `V` is a tuple or vector of variables and `P` a tuple or vector of polynomials, e.g. `(x, y) => (y, x)` or `(y, x) => (2, 1)`.

The order of the variables is lexicographic with the name with TypedPolynomials and by order of creation with DynamicPolynomials.
Since there is no guarantee on the order of the variables, substitution directly with a tuple or a vetor is not allowed.
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

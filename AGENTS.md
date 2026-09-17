
This repository as well as DynamicPolynomials have accumulated a lot of technical depth due to bad design decisions and it is time to get rid of it.
The bad design decisions were
1. Try to support `*(::MP.AbstractPolynomialLike, ::Any)` this creates a lot of invalidation and promotion rules are a mess. We should only support `*(::MP.AbstractPolynomialLike{T}, ::T)` and have `*(::MP.AbstractPolynomialLike, ::Number)` that convert the number to `T`. Same for other operators like `+`, `-` etc..
2. Implement `+(p::MP.AbstractPolynomialLike, q::MP.AbstractPolynomialLike)` where `p` and `q` have different variables. Not simply by first promoting them to the same set of variables and then use a simple implementation of the sum but instead to have some complicated implementation dealing with monomials over different variables. We should just directly promote with `StarAlgebras.promote_bases` now and only implement `+` over polynomials of the same algebra
3. Having different implementation of polynomials, `MP.Polynomial` using a list of terms and `DynamicPolynomials.Polynomial` having a separate list of coefficients and list of monomials. We should just use StarAlgebras.AlgebraElement now.

So the plan is:
1. Move `MP.Term` to `StarAlgebras.Term`, so `MP.AbstractTermLike` won't be an abstract type anymore but a union of `AbstractMonomialLike` and `SA.Term`
2. Remove `MP.Polynomial`, and replace it by the implementation of polynomials done in MultivariateBases using StarAlgebras. Doing so, we can just merge MultivariateBases into MultivariatePolynomials, `MP.AbstractPolynomialLike` will then be a union of that StarAlgebras.AlgebraElement and `AbstractTermLike`.

Terminology:

MOI: MathOptInterface
MA: MutableArithmetics
DP: DynamicPolynomials
TP: TypedPolynomials
MP: MultivariatePolynomials
SA: StarAlgebras

All packages should be in ~/.julia/dev

MA has a large test suite that is used by JuMP, MOI, MP and therefore also DP and TP. It is in MP/test. It's nice but it takes a lot of time to run. First makes sure everything is done before running it for the final touches.
This refactor is supposed to get rid of a lot of code and technical depth, don't start writing a lot of code before first asking me.

The promotion rules are a mess. Because we had to support promotion with Any, this messed up badly with Julia's internals.
In Julia, normally you just need to implement promote_rule(::A, ::B), not promote_rule(::B, ::A). But because I implemented promotion with ::Any, I had to do both and it created such a mess!!
It was also a rabbithole where I then needed also APL{<:Any} and things like that, quite complicated...
I want things to be easier now. Before, I belived if x is a variable, promote_type(typeof(x), typeof("x")) would be something like Term{String}. We don't want that anymore, let it just be Any[].
You see the vibe ? Let's start simple and then be really careful for test we used to have whether we really want this test to still pass or we just want to simplify.

About backward compat, we are going to make a breaking release so I prefer being more breaking and having a simpler code than the opposite.

A lot of code in MP are like sum(t::Vector{<:AbstractTerm}) = dot(coefficient.(t), monomial.(t)) and then dot(c::AbstractVector, m::AbstractVector{<:AbstractMonomialLike}) = sum(c .* m).
These create stackoverflow in case DP or TP forget to implement one of the two! This was because, as explained above, TP used the default polynomials that were a vector of terms and DP
was a vector of coefficients separated to a vector of monomials. Now it all gets simpler because it will be a separate vector of coefficients and monomials since we'll just be using SA.AlgebraElement!!
So we can simplify. Also, because SA.AlgebraElement <: MA.AbstractMutable, we already have a fallback for Base.sum so we might just be able to remove these!
These were written before MutableArithmetics existed. I am maintaining MA, SA, JuMP, MOI, MP, DP, TP and I don't want to maintain duplicates. All these packages define mutable objects.
MA allows me to have already a lot of code in common of all of them. Then, JuMP, SA and MOI will have different implementations of basically "sum of terms" but that's fine.
What I don't want is MP, TP and DP to also have their own version, they should just use SA.AlgebraElement!

It has already been started by another AI agent on the branches:
DP: bl/sa_term
MP: bl/sa_poly
SA: bl/term

Just continue on these branches, feel free to simplify them if you see anything better, it was done a while ago by an older AI, you are smarter ;)

In the printing of stack-traces, you can see that because we just have generic types in SA that we parametrize, things are getting very long and it's getting very difficult to debug.
This gives the user extra flexibility to try new these, and in these cases, it will be nice to have precise stacktraces, but we also want the common cases (like what you get with DP.@polyvar x y; 2 * x + y)
to have a very small type like DP.Polynomial{Int}
We can solve it easily by having a "const Polynomial = ..." in DynamicPolynomials, don't hesitate to do this early on, it will help you be more context-efficient.

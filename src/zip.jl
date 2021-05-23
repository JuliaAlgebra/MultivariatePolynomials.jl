"""
pairzip((a, b), (c, d)) gives (a=>c, b=>d)

This function was written by Fengyang Wang and shared on the Julia discourse
forum: https://discourse.julialang.org/t/type-stable-zip-to-pairs/3390/2
"""
pairzip(::Tuple{}, ::Tuple{}) = ()
pairzip(::Tuple{}, ::Tuple) = throw(ArgumentError("args must be equal in length"))
pairzip(::Tuple, ::Tuple{}) = throw(ArgumentError("args must be equal in length"))
pairzip(t::Tuple, u::Tuple) = (t[1] => u[1], pairzip(Base.tail(t), Base.tail(u))...)
pairzip(p::Pair{<:Tuple, <:Tuple}) = pairzip(p.first, p.second)
# inefficient but convenient method to allow subs(p, (x, y)=>[1, 2])
pairzip(p::Pair) = pairzip(Tuple(p.first), Tuple(p.second))

tuplezip(::Tuple{}, ::Tuple{}) = ()
tuplezip(::Tuple{}, ::Tuple) = throw(ArgumentError("args must be equal in length"))
tuplezip(::Tuple, ::Tuple{}) = throw(ArgumentError("args must be equal in length"))
tuplezip(t::Tuple, u::Tuple) = ((t[1], u[1]), tuplezip(Base.tail(t), Base.tail(u))...)

_zip(t::Tuple, u::Tuple) = tuplezip(t, u)
# `tuplezip` would be type unstable
_zip(t::Vector, u::Vector) = zip(t, u)

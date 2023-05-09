"""
pair_zip((a, b), (c, d)) gives (a=>c, b=>d)

This function was written by Fengyang Wang and shared on the Julia discourse
forum: https://discourse.julialang.org/t/type-stable-zip-to-pairs/3390/2
"""
pair_zip(::Tuple{}, ::Tuple{}) = ()
pair_zip(::Tuple{}, ::Tuple) = throw(ArgumentError("args must be equal in length"))
pair_zip(::Tuple, ::Tuple{}) = throw(ArgumentError("args must be equal in length"))
pair_zip(t::Tuple, u::Tuple) = (t[1] => u[1], pair_zip(Base.tail(t), Base.tail(u))...)
pair_zip(p::Pair{<:Tuple, <:Tuple}) = pair_zip(p.first, p.second)
# inefficient but convenient method to allow subs(p, (x, y)=>[1, 2])
pair_zip(p::Pair) = pair_zip(Tuple(p.first), Tuple(p.second))

tuple_zip(::Tuple{}, ::Tuple{}) = ()
tuple_zip(::Tuple{}, ::Tuple) = throw(ArgumentError("args must be equal in length"))
tuple_zip(::Tuple, ::Tuple{}) = throw(ArgumentError("args must be equal in length"))
tuple_zip(t::Tuple, u::Tuple) = ((t[1], u[1]), tuple_zip(Base.tail(t), Base.tail(u))...)

_zip(t::Tuple, u::Tuple) = tuple_zip(t, u)
# `tuple_zip` would be type unstable
_zip(t::Vector, u::Vector) = zip(t, u)

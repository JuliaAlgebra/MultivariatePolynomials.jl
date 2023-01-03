# Allows `terms(::AbstractTerm)` to not allocate and be type stable
struct OneOrZeroElementVector{T} <: AbstractVector{T}
    empty::Bool
    el::T
end

Base.size(it::OneOrZeroElementVector) = (length(it),)

Base.length(it::OneOrZeroElementVector) = Int(!it.empty)

function Base.iterate(it::OneOrZeroElementVector)
    if it.empty
        return nothing
    else
        return it.el, nothing
    end
end
Base.iterate(::OneOrZeroElementVector, ::Nothing) = nothing

Base.getindex(it::OneOrZeroElementVector, i::Integer) = it.el

Iterators.reverse(it::OneOrZeroElementVector) = it

# Copied from MathOptInterface
"""
    struct LazyMap{T, VT}
        f::Function
        data::VT
    end

Iterator over the elements of `data` mapped by `f`. This is similar to
`Base.Generator(f, data)` except that the `eltype` of a `LazyMap` is given at
construction while the `eltype` of `Base.Generator(f, data)` is `Any`.
"""
struct LazyMap{T,VT,F}
    f::F
    data::VT
end

function LazyMap{T}(f, data) where {T}
    return LazyMap{T,typeof(data),typeof(f)}(f, data)
end

Base.size(it::LazyMap) = size(it.data)

Base.length(it::LazyMap) = length(it.data)

function Base.iterate(it::LazyMap, args...)
    elem_state = iterate(it.data, args...)
    if elem_state === nothing
        return
    else
        return it.f(elem_state[1]), elem_state[2]
    end
end

Base.IteratorSize(it::LazyMap) = Base.IteratorSize(it.data)

Base.eltype(::LazyMap{T}) where {T} = T

Base.getindex(it::LazyMap, i) = it.f(getindex(it.data, i))

Base.eachindex(it::LazyMap) = Base.eachindex(it.data)
Base.lastindex(it::LazyMap) = Base.lastindex(it.data)

function Iterators.reverse(it::LazyMap{T}) where {T}
    return LazyMap{T}(it.f, Iterators.reverse(it.data))
end

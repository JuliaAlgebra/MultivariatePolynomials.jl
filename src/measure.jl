export Measure, zeta, ζ

struct Moment{T, MT <: AbstractMonomial}
    α::T
    x::MT
end

# If a monomial is not in x, it does not mean that the moment is zero, it means that it is unknown/undefined
struct Measure{T, MT <: AbstractMonomial, MVT <: AbstractVector{MT}}
    a::Vector{T}
    x::MVT

    function Measure{T, MT, MVT}(a::Vector{T}, x::MVT) where {T, MT, MVT}
        @assert length(a) == length(x)
        new(a, x)
    end
end

Measure(a::Vector{T}, x::AbstractVector{MT}) where {T, MT <: AbstractMonomial} = Measure{T, MT, monovectype(x)}(monovec(a, x)...)

function ζ(x::AbstractVector{MT}, s::AbstractSubstitution...) where {MT <: AbstractMonomial}
    Measure([m(s...) for m in x], x)
end

struct MatMeasure{T, MT <: AbstractMonomial, MVT <: AbstractVector{MT}}
    Q::Vector{T}
    x::MVT
end

export Measure, zeta, ζ

type Moment{T, MT <: AbstractMonomial}
    α::T
    x::MT
end

# If a monomial is not in x, it does not mean that the moment is zero, it means that it is unknown/undefined
type Measure{T, MT <: AbstractMonomial, MVT <: AbstractVector{MT}}
    a::Vector{T}
    x::MVT

    function Measure{T, MT, MVT}(a::Vector{T}, x::MVT) where {T, MT, MVT}
        @assert length(a) == length(x)
        new(a, x)
    end
end

Measure{T, MT <: AbstractMonomial}(a::Vector{T}, x::AbstractVector{MT}) = Measure{T, MT, monovectype(x)}(monovec(a, x)...)

function ζ{T, MT <: AbstractMonomial, PVT <: AbstractVariable}(v::Vector{T}, x::AbstractVector{MT}, varorder::AbstractVector{PVT})
    Measure(T[m(v, varorder) for m in x], x)
end

type MatMeasure{T, MT <: AbstractMonomial, MVT <: AbstractVector{MT}}
    Q::Vector{T}
    x::MVT
end

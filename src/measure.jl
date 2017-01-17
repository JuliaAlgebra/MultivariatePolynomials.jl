export Measure, zeta, ζ

abstract AbstractMeasure

# If a monomial is not in x, it does not mean that the moment is zero, it means that it is unknown/undefined
type Measure{C, T}
    a::Vector{T}
    x::MonomialVector{C}

    function Measure(a::Vector{T}, x::MonomialVector{C})
        if length(a) != length(x)
            error("There should be as many coefficient than monomials")
        end
        new(a, x)
    end
end

Measure{C, T}(a::Vector{T}, x::MonomialVector{C}) = Measure{C, T}(a, x)
function (::Type{Measure{C}}){C}(a::Vector, x::Vector)
    if length(a) != length(x)
        error("There should be as many coefficient than monomials")
    end
    perm, X = sortmonovec(PolyVar{C}, x)
    Measure(a[perm], X)
end
Measure{T<:VectorOfPolyType{true}}(a::Vector, X::Vector{T}) = Measure{true}(a, X)
Measure{T<:VectorOfPolyType{false}}(a::Vector, X::Vector{T}) = Measure{false}(a, X)

function ζ{C, T}(v::Vector{T}, x::MonomialVector{C}, varorder::Vector{PolyVar{C}})
    Measure(T[m(v, varorder) for m in x], x)
end

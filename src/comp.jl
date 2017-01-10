import Base.==, Base.isless, Base.isapprox

# TODO this should be in Base
function (==){T}(x::Vector{T}, y::Vector{T})
    if length(x) != length(y)
        false
    else
        #for (xi, yi) in zip(x, y)
        for i in 1:length(x)
            if x[i] != y[i]
                return false
            end
        end
        true
    end
end

# Technique: the higher catch the calls when it is rhs
# so p::PolyType == x::PolyVar -> x == p
(==)(p::PolyType, y) = y == p

# TODO equality should be between name ?
function (==){PV<:AbstractPolyVar}(x::PV, y::PV)
    x === y
end
(==)(α, x::AbstractPolyVar) = false
(==)(p::PolyType, x::AbstractPolyVar) = x == p

isless{PV<:AbstractPolyVar}(x::PV, y::PV) = isless(y.name, x.name)

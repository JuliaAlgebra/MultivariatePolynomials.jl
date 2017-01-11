export NCPolyVar, NCMonomial, NCMonomialVector, @ncpolyvar

macro ncpolyvar(args...)
    reduce((x,y) -> :($x; $y), :(), [buildpolyvar(arg, false) for arg in args])
end

immutable NCPolyVar <: PolyType
    name::AbstractString
end

type NCMonomial <: PolyType
    vars::Vector{NCPolyVar}
    z::Vector{Int}

    function NCMonomial(vars::Vector{NCPolyVar}, z::Vector{Int})
        if length(vars) != length(z)
            throw(ArgumentError("There should be as many vars than exponents"))
        end
        new(vars, z)
    end
end
NCMonomial() = NCMonomial(NCPolyVar[], Int[])
Base.convert(::Type{NCMonomial}, x::NCPolyVar) = NCMonomial([x], [1])

type NCMonomialVector <: PolyType
    vars::Vector{NCPolyVar}
    Z::Vector{Vector{Int}}

    function NCMonomialVector(vars::Vector{NCPolyVar}, Z::Vector{Vector{Int}})
        for z in Z
            if length(vars) != length(z)
                throw(ArgumentError("There should be as many vars than exponents"))
            end
        end
        @assert issorted(Z, rev=true)
        new(vars, Z)
    end
end
NCMonomialVector() = NCMonomialVector(NCPolyVar[], Vector{Int}[])

function getindex(x::NCMonomialVector, i::Integer)
    NCMonomial(x.vars, x.Z[i])
end

function sortmonovec{T<:Union{PolyType,Int}}(::Type{NCPolyVar}, X::Vector{T})
    allvars, Z = buildZvarsvec(NCPolyVar, X)
    perm = sortperm(Z, rev=true)
    perm, NCMonomialVector(allvars, Z[perm])
end
function NCMonomialVector{T<:Union{PolyType,Int}}(X::Vector{T})
    allvars, Z = buildZvarsvec(NCPolyVar, X)
    sort!(Z, rev=true)
    NCMonomialVector(allvars, Z)
end

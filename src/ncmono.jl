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

type NCMonomialVector <: PolyType
    vars::Vector{NCPolyVar}
    Z::Vector{Vector{Int}}

    function NCMonomialVector(vars::Vector{NCPolyVar}, Z::Vector{Vector{Int}})
        for z in Z
            if length(vars) != length(z)
                error("There should be as many vars than exponents")
            end
        end
        @assert issorted(Z, rev=true)
        new(vars, Z)
    end
end
NCMonomialVector() = NCMonomialVector(NCPolyVar[], Vector{Int}[])

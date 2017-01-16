export NCPolyVar, NCMonomial, NCMonomialVector, @ncpolyvar

macro ncpolyvar(args...)
    reduce((x,y) -> :($x; $y), :(), [buildpolyvar(arg, false) for arg in args])
end

immutable NCPolyVar <: AbstractPolyVar
    name::AbstractString
end

type NCMonomial <: AbstractMonomial
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

type NCMonomialVector <: AbstractMonomialVector
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

function fillZrec!(Z, z, i, n, deg, filter::Function)
    if deg == 0
        if filter(z)
            push!(Z, copy(z))
        end
    else
        for i in i:i+n-1
            z[i] += 1
            fillZrec!(Z, z, i, n, deg-1, filter)
            z[i] -= 1
        end
    end
end
function fillZfordeg!(Z, n, deg, ::Type{Val{false}}, filter::Function)
    z = zeros(Int, deg * n - deg + 1)
    fillZrec!(Z, z, 1, n, deg, filter)
end

function getvarsforlength(vars::Vector{NCPolyVar}, len::Int)
    n = length(vars)
    map(i -> vars[((i-1) % n) + 1], 1:len)
end
function NCMonomialVector(vars::Vector{NCPolyVar}, degs, filter::Function = x->true)
    Z = getZfordegs(length(vars), degs, false, filter)
    v = isempty(Z) ? vars : getvarsforlength(vars, length(first(Z)))
    NCMonomialVector(v, Z)
end
NCMonomialVector(vars::Vector{NCPolyVar}, degs::Int, filter::Function = x->true) = NCMonomialVector(vars, [degs], filter)
function monomials(vars::Vector{NCPolyVar}, degs::Vector{Int}, filter::Function = x->true)
    Z = getZfordegs(length(vars), degs, false, filter)
    v = isempty(Z) ? vars : getvarsforlength(vars, length(first(Z)))
    [NCMonomial(v, z) for z in Z]
end

function sortmonovec{T<:Union{NCMonomial,NCPolyVar,Int}}(::Type{NCPolyVar}, X::Vector{T})
    allvars, Z = buildZvarsvec(NCPolyVar, X)
    perm = sortperm(Z, rev=true)
    perm, NCMonomialVector(allvars, Z[perm])
end
function NCMonomialVector{T<:Union{NCMonomial,NCPolyVar,Int}}(X::Vector{T})
    allvars, Z = buildZvarsvec(NCPolyVar, X)
    sort!(Z, rev=true)
    NCMonomialVector(allvars, Z)
end

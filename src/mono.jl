export PolyVar, Monomial, MonomialVector, @polyvar, @ncpolyvar
export monomials, polyvecvar, vars, nvars, extdeg, mindeg, maxdeg

function polyvecvar(prefix, idxset, C::Bool)
    [PolyVar{C}("$(prefix * string(i))") for i in idxset]
end

function buildpolyvar(var, C::Bool)
    if isa(var, Symbol)
        :($(esc(var)) = PolyVar{$C}($"$var"))
    else
        isa(var, Expr) || error("Expected $var to be a variable name")
        Base.Meta.isexpr(var, :ref) || error("Expected $var to be of the form varname[idxset]")
        length(var.args) == 2 || error("Expected $var to have one index set")
        varname = var.args[1]
        prefix = string(var.args[1])
        idxset = esc(var.args[2])
        :($(esc(varname)) = polyvecvar($prefix, $idxset, true))
    end
end

# Variable vector x returned garanteed to be sorted so that if p is built with x then vars(p) == x
macro polyvar(args...)
    reduce((x,y) -> :($x; $y), :(), [buildpolyvar(arg, true) for arg in args])
end
macro ncpolyvar(args...)
    reduce((x,y) -> :($x; $y), :(), [buildpolyvar(arg, false) for arg in args])
end

# TODO equality should be between name ?
immutable PolyVar{C} <: PolyType{C}
    name::AbstractString
end
iscomm{C}(::Type{PolyVar{C}}) = C

copy(x::PolyVar) = x

vars(x::PolyVar) = [x]
nvars(::PolyVar) = 1

function myunion{PV<:PolyVar}(varsvec::Vector{Vector{PV}})
    n = length(varsvec)
    is = ones(Int, n)
    maps = [ zeros(Int, length(vars)) for vars in varsvec ]
    nonempty = IntSet(find([!isempty(vars) for vars in varsvec]))
    vars = Vector{PV}()
    while !isempty(nonempty)
        imin = 0
        for i in nonempty
            if imin == 0 || varsvec[i][is[i]] > varsvec[imin][is[imin]]
                imin = i
            end
        end
        var = varsvec[imin][is[imin]]
        push!(vars, var)
        for i in nonempty
            if var == varsvec[i][is[i]]
                maps[i][is[i]] = length(vars)
                if is[i] == length(varsvec[i])
                    pop!(nonempty, i)
                else
                    is[i] += 1
                end
            end
        end
    end
    vars, maps
end

# Invariant:
# vars is increasing
# z may contain 0's (otherwise, getindex of MonomialVector would be inefficient)
type Monomial{C} <: PolyType{C}
    vars::Vector{PolyVar{C}}
    z::Vector{Int}

    function Monomial(vars::Vector{PolyVar{C}}, z::Vector{Int})
        if length(vars) != length(z)
            throw(ArgumentError("There should be as many vars than exponents"))
        end
        new(vars, z)
    end
end
iscomm{C}(::Type{Monomial{C}}) = C
(::Type{Monomial{C}}){C}() = Monomial{C}(PolyVar{C}[], Int[])
Base.convert{C}(::Type{Monomial{C}}, x::PolyVar{C}) = Monomial{C}([x], [1])
Monomial{C}(vars::Vector{PolyVar{C}}, z::Vector{Int}) = Monomial{C}(vars, z)
Monomial{C}(x::PolyVar{C}) = Monomial{C}(x)

# /!\ vars not copied, do not mess with vars
copy{M<:Monomial}(m::M) = M(m.vars, copy(m.z))
deg(x::Monomial) = sum(x.z)
nvars(x::Monomial) = length(x.vars)
isconstant(x::Monomial) = deg(x) == 0

# Invariant: Always sorted and no zero vector
type MonomialVector{C} <: PolyType{C}
    vars::Vector{PolyVar{C}}
    Z::Vector{Vector{Int}}

    function MonomialVector(vars::Vector{PolyVar{C}}, Z::Vector{Vector{Int}})
        for z in Z
            if length(vars) != length(z)
                throw(ArgumentError("There should be as many vars than exponents"))
            end
        end
        @assert issorted(Z, rev=true)
        new(vars, Z)
    end
end
MonomialVector{C}(vars::Vector{PolyVar{C}}, Z::Vector{Vector{Int}}) = MonomialVector{C}(vars, Z)
(::Type{MonomialVector{C}}){C}() = MonomialVector{C}(PolyVar{C}[], Vector{Int}[])

# /!\ vars not copied, do not mess with vars
copy{MV<:MonomialVector}(m::MV) = MV(m.vars, copy(m.Z))
function getindex{MV<:MonomialVector}(x::MV, I)
    MV(x.vars, x.Z[I])
end

length(x::MonomialVector) = length(x.Z)
isempty(x::MonomialVector) = length(x) == 0
start(::MonomialVector) = 1
done(x::MonomialVector, state) = length(x) < state
next(x::MonomialVector, state) = (x[state], state+1)

extdeg(x::MonomialVector) = extrema(sum.(x.Z))
mindeg(x::MonomialVector) = minimum(sum.(x.Z))
maxdeg(x::MonomialVector) = maximum(sum.(x.Z))

vars{T<:Union{Monomial, MonomialVector}}(x::T) = x.vars


function getindex(x::MonomialVector, i::Integer)
    Monomial(x.vars, x.Z[i])
end

function fillZfordeg!(Z, n, deg, ::Type{Val{true}}, filter::Function)
    z = zeros(Int, n)
    z[1] = deg
    while true
        if filter(z)
            push!(Z, z)
            z = copy(z)
        end
        if z[end] == deg
            break
        end
        sum = 1
        for j in (n-1):-1:1
            if z[j] != 0
                z[j] -= 1
                z[j+1] += sum
                break
            else
                sum += z[j+1]
                z[j+1] = 0
            end
        end
    end
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
# List exponents in decreasing Graded Lexicographic Order
function getZfordegs(n, degs::AbstractVector{Int}, C::Bool, filter::Function)
    Z = Vector{Vector{Int}}()
    for deg in sort(degs, rev=true)
        fillZfordeg!(Z, n, deg, Val{C}, filter)
    end
    @assert issorted(Z, rev=true)
    Z
end
function MonomialVector(vars::Vector{PolyVar{true}}, degs::AbstractVector{Int}, filter::Function = x->true)
    MonomialVector{true}(vars, getZfordegs(length(vars), degs, true, filter))
end
function getvarsforlength(vars::Vector{PolyVar{false}}, len::Int)
    n = length(vars)
    map(i -> vars[((i-1) % n) + 1], 1:len)
end

function MonomialVector(vars::Vector{PolyVar{false}}, degs::AbstractVector{Int}, filter::Function = x->true)
    Z = getZfordegs(length(vars), degs, false, filter)
    v = isempty(Z) ? vars : getvarsforlength(vars, length(first(Z)))
    MonomialVector{false}(v, Z)
end
MonomialVector{C}(vars::Vector{PolyVar{C}}, degs::Int, filter::Function = x->true) = MonomialVector(vars, [degs], filter)
function monomials(vars::Vector{PolyVar{true}}, degs::AbstractVector{Int}, filter::Function = x->true)
    Z = getZfordegs(length(vars), degs, true, filter)
    [Monomial{true}(vars, z) for z in Z]
end
function monomials(vars::Vector{PolyVar{false}}, degs::AbstractVector{Int}, filter::Function = x->true)
    Z = getZfordegs(length(vars), degs, false, filter)
    v = isempty(Z) ? vars : getvarsforlength(vars, length(first(Z)))
    [Monomial{false}(v, z) for z in Z]
end
monomials{PV<:PolyVar}(vars::Vector{PV}, degs::Int, filter::Function = x->true) = monomials(vars, [degs], filter)

function buildZvarsvec{PV<:PolyVar, T<:Union{PolyType,Int}}(::Type{PV}, X::Vector{T})
    varsvec = Vector{PV}[ (isa(x, PolyType) ? vars(x) : PolyVar[]) for x in X ]
    allvars, maps = myunion(varsvec)
    nvars = length(allvars)
    Z = [zeros(Int, nvars) for i in 1:length(X)]
    offset = 0
    for (i, x) in enumerate(X)
        if isa(x, PolyVar)
            @assert length(maps[i]) == 1
            z = [1]
        elseif isa(x, Monomial)
            z = x.z
        elseif isa(x, Term)
            z = x.x.z
        else
            @assert isa(x, Int)
            z = Int[]
        end
        Z[i][maps[i]] = z
    end
    allvars, Z
end
function sortmonovec{C, T<:Union{PolyType,Int}}(::Type{PolyVar{C}}, X::Vector{T})
    allvars, Z = buildZvarsvec(PolyVar{C}, X)
    perm = sortperm(Z, rev=true)
    perm, MonomialVector{C}(allvars, Z[perm])
end
function (::Type{MonomialVector{C}}){C}(X::Vector)
    allvars, Z = buildZvarsvec(PolyVar{C}, X)
    sort!(Z, rev=true)
    MonomialVector{C}(allvars, Z)
end
typealias VectorOfPolyType{C} Union{PolyType{C},Int}
MonomialVector{T<:VectorOfPolyType{false}}(X::Vector{T}) = MonomialVector{false}(X)
MonomialVector{T<:VectorOfPolyType{true}}(X::Vector{T}) = MonomialVector{true}(X)

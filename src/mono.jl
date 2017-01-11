export monomials, polyvecvar, vars, nvars, extdeg, mindeg, maxdeg

function polyvecvar(prefix, idxset, comm)
    PV = comm ? PolyVar : NCPolyVar
    [PV("$(prefix * string(i))") for i in idxset]
end

function buildpolyvar(var, comm)
    if isa(var, Symbol)
        PV = comm ? PolyVar : NCPolyVar
        :($(esc(var)) = $PV($"$var"))
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

typealias AbstractPolyVar Union{PolyVar, NCPolyVar}
copy(x::AbstractPolyVar) = x

vars(x::AbstractPolyVar) = [x]
nvars(::AbstractPolyVar) = 1

function myunion{PV<:AbstractPolyVar}(varsvec::Vector{Vector{PV}})
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

typealias AbstractMonomial Union{Monomial, NCMonomial}
# /!\ vars not copied, do not mess with vars
copy{M<:AbstractMonomial}(m::M) = M(m.vars, copy(m.z))
deg(x::AbstractMonomial) = sum(x.z)
nvars(x::AbstractMonomial) = length(x.vars)
isconstant(x::AbstractMonomial) = deg(x) == 0

typealias AbstractMonomialVector Union{MonomialVector, NCMonomialVector}
# /!\ vars not copied, do not mess with vars
copy{MV<:AbstractMonomialVector}(m::MV) = MV(m.vars, copy(m.Z))
function getindex{MV<:AbstractMonomialVector}(x::MV, I)
    MV(x.vars, x.Z[I])
end

length(x::AbstractMonomialVector) = length(x.Z)
isempty(x::AbstractMonomialVector) = length(x) == 0
start(::AbstractMonomialVector) = 1
done(x::AbstractMonomialVector, state) = length(x) < state
next(x::AbstractMonomialVector, state) = (x[state], state+1)

extdeg(x::AbstractMonomialVector) = extrema(sum.(x.Z))
mindeg(x::AbstractMonomialVector) = minimum(sum.(x.Z))
maxdeg(x::AbstractMonomialVector) = maximum(sum.(x.Z))

vars{T<:Union{AbstractMonomial, AbstractMonomialVector}}(x::T) = x.vars

function buildZvarsvec{PV<:AbstractPolyVar, T<:Union{PolyType,Int}}(::Type{PV}, X::Vector{T})
    varsvec = Vector{PV}[ (isa(x, PolyType) ? vars(x) : PolyVar[]) for x in X ]
    allvars, maps = myunion(varsvec)
    nvars = length(allvars)
    Z = [zeros(Int, nvars) for i in 1:length(X)]
    offset = 0
    for (i, x) in enumerate(X)
        if isa(x, AbstractPolyVar)
            @assert length(maps[i]) == 1
            z = [1]
        elseif isa(x, AbstractMonomial)
            z = x.z
        elseif isa(x, AbstractTerm)
            z = x.x.z
        else
            @assert isa(x, Int)
            z = Int[]
        end
        Z[i][maps[i]] = z
    end
    allvars, Z
end

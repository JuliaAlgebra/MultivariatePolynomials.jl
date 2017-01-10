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

^(x::NCPolyVar, i::Int) = NCMonomial([x], [i])

function (*)(x::NCPolyVar, y::NCPolyVar)
    if x === y
        NCMonomial([x], [2])
    else
        NCMonomial(x > y ? [x,y] : [y,x], [1,1])
    end
end

function multiplyvar(v::Vector{NCPolyVar}, x::NCPolyVar)
    if v[end] === x
        multiplyexistingvar(v, x, length(v))
    else
        insertvar(v, x, length(v)+1)
    end
end
function multiplyvar(x::NCPolyVar, v::Vector{NCPolyVar})
    if v[1] === x
        multiplyexistingvar(v, x, 1)
    else
        insertvar(v, x, 1)
    end
end

function (*)(x::NCPolyVar, y::NCMonomial)
    w, updatez = multiplyvar(x, y.vars)
    NCMonomial(w, updatez(y.z))
end
function (*)(y::NCMonomial, x::NCPolyVar)
    w, updatez = multiplyvar(y.vars, x)
    NCMonomial(w, updatez(y.z))
end
function (*)(x::NCPolyVar, y::NCMonomialVector)
    w, updatez = multiplyvar(x, y.vars)
    NCMonomialVector(w, updatez.(y.Z))
end
function (*)(y::NCMonomialVector, x::NCPolyVar)
    w, updatez = multiplyvar(y.vars, x)
    NCMonomialVector(w, updatez.(y.Z))
end
function (*)(x::NCMonomial, y::NCMonomial)
    w, updatez = multiplymono(y.vars, x)
    NCMonomial(w, updatez(y.z))
end
function (*)(x::NCMonomial, y::NCMonomialVector)
    w, updatez = multiplymono(y.vars, x)
    NCMonomialVector(w, updatez.(y.Z))
end
function (*)(y::NCMonomialVector, x::NCMonomial)
    w, updatez = multiplymono(y.vars, x)
    NCMonomialVector(w, updatez.(y.Z))
end

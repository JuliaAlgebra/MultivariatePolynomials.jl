^(x::NCPolyVar, i::Int) = NCMonomial([x], [i])

function (*)(x::NCPolyVar, y::NCPolyVar)
    if x === y
        NCMonomial([x], [2])
    else
        NCMonomial([x, y], [1, 1])
    end
end

function multiplyvar(v::Vector{NCPolyVar}, z::Vector{Int}, x::NCPolyVar)
    i = length(v)
    while i > 0 && z[i] == 0
        i -= 1
    end
    if v[i] == x
        multiplyexistingvar(v, x, i)
    else
        #   ---->
        # \  |\  |\
        #  \ | \ | \
        #   \|  \|  \
        # If z[i] > x, we wait either for a rise (v[i] > v[i-1]) or v[i] < x
        # Otherwise, we first wait for a drop and then wait for the same thing
        ndrop = 0
        if v[i] > x
            droplim1 = 0
            droplim2 = 1
        else
            droplim1 = 1
            droplim2 = 2
        end
        i += 1
        while i <= length(v) && v[i] != x
            if v[i] > v[i-1]
                ndrop += 1
            end
            if drop >= droplim2 || (drop >= droplim1 && v[i] < x)
                break
            end
            i += 1
        end

        if i <= length(v) && v[i] == x
            multiplyexistingvar(v, x, i)
        else
            insertvar(v, x, i)
        end
    end
end
function multiplyvar(x::NCPolyVar, v::Vector{NCPolyVar}, z::Vector{Int})
    i = 1
    while i <= length(v) && z[i] == 0
        i += 1
    end
    if v[i] == x
        multiplyexistingvar(v, x, i)
    else
        #   <----
        # \  |\  |\
        #  \ | \ | \
        #   \|  \|  \
        # If z[i] < x, we wait either for a drop (v[i] < v[i+1]) or v[i] > x
        # Otherwise, we first wait for a drop and then wait for the same thing
        ndrop = 0
        if v[i] < x
            droplim1 = 0
            droplim2 = 1
        else
            droplim1 = 1
            droplim2 = 2
        end
        i -= 1
        while i > 0 && v[i] != x
            if v[i] < v[i+1]
                ndrop += 1
            end
            if drop >= droplim2 || (drop >= droplim1 && v[i] > x)
                break
            end
            i -= 1
        end
        if i > 0 && v[i] == x
            multiplyexistingvar(v, x, i)
        else
            insertvar(v, x, i+1)
        end
    end
end
function (*)(x::NCPolyVar, y::NCMonomial)
    w, updatez = multiplyvar(x, y.vars, y.z)
    NCMonomial(w, updatez(y.z))
end
function (*)(y::NCMonomial, x::NCPolyVar)
    w, updatez = multiplyvar(y.vars, y.z, x)
    NCMonomial(w, updatez(y.z))
end

function multiplymono(v::Vector{NCPolyVar}, x::NCMonomial)
    w, maps = myunion([v, x.vars])
    updatez = z -> begin
        newz = zeros(Int, length(w))
        newz[maps[1]] += z
        newz[maps[2]] += x.z
        newz
    end
    w, updatez
end
function (*)(x::NCMonomial, y::NCMonomial)
    i = findlast(z -> z > 0, x.z)
    if i == 0
        return y
    end
    j = findfirst(z -> z > 0, y.z)
    if j == 0
        return x
    end
    if x.vars[i] == y.vars[j]
        w = [x.vars[1:i]; y.vars[j+1:end]]
        z = [x.z[1:i-1]; x.z[i] + y.z[j]; y.z[j+1:end]]
    else
        w = [x.vars[1:i]; y.vars[j:end]]
        z = [x.z[1:i]; y.z[j:end]]
    end
    return NCMonomial(w, z)
end

function (*)(x::NCPolyVar, y::NCMonomialVector)
    NCMonomialVector([x * yi for yi in y])
end
function (*)(y::NCMonomialVector, x::NCPolyVar)
    NCMonomialVector([yi * x for yi in y])
end
function (*)(x::NCMonomial, y::NCMonomialVector)
    NCMonomialVector([x * yi for yi in y])
end
function (*)(y::NCMonomialVector, x::NCMonomial)
    NCMonomialVector([yi * x for yi in y])
end

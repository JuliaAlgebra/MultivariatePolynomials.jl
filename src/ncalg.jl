^(x::PolyVar{false}, i::Int) = Monomial{false}([x], [i])

function (*)(x::PolyVar{false}, y::PolyVar{false})
    if x === y
        Monomial{false}([x], [2])
    else
        Monomial{false}([x, y], [1, 1])
    end
end

function multiplyvar(v::Vector{PolyVar{false}}, z::Vector{Int}, x::PolyVar{false})
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
            if ndrop >= droplim2 || (ndrop >= droplim1 && v[i] < x)
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
function multiplyvar(x::PolyVar{false}, v::Vector{PolyVar{false}}, z::Vector{Int})
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
            if ndrop >= droplim2 || (ndrop >= droplim1 && v[i] > x)
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
function (*)(x::PolyVar{false}, y::Monomial{false})
    w, updatez = multiplyvar(x, y.vars, y.z)
    Monomial{false}(w, updatez(y.z))
end
function (*)(y::Monomial{false}, x::PolyVar{false})
    w, updatez = multiplyvar(y.vars, y.z, x)
    Monomial{false}(w, updatez(y.z))
end

function (*)(x::Monomial{false}, y::Monomial{false})
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
    return Monomial{false}(w, z)
end

function (*)(x::PolyVar{false}, y::MonomialVector{false})
    MonomialVector{false}([x * yi for yi in y])
end
function (*)(y::MonomialVector{false}, x::PolyVar{false})
    MonomialVector{false}([yi * x for yi in y])
end
function (*)(x::Monomial{false}, y::MonomialVector{false})
    MonomialVector{false}([x * yi for yi in y])
end
function (*)(y::MonomialVector{false}, x::Monomial{false})
    MonomialVector{false}([yi * x for yi in y])
end

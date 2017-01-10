function multiplyexistingvar{PV<:AbstractPolyVar}(v::Vector{PV}, x::PV, i::Int)
    updatez = z -> begin
        newz = copy(z)
        newz[i] += 1
        newz
    end
    # /!\ v not copied for efficiency, do not mess up with vars
    v, updatez
end
function insertvar{PV<:AbstractPolyVar}(v::Vector{PV}, x::PV, i::Int)
    n = length(v)
    I = 1:i-1
    J = i:n
    K = J+1
    w = Vector{PV}(n+1)
    w[I] = v[I]
    w[i] = x
    w[K] = v[J]
    updatez = z -> begin
        newz = Vector{Int}(n+1)
        newz[I] = z[I]
        newz[i] = 1
        newz[K] = z[J]
        newz
    end
    w, updatez
end

function multiplymono{PV<:AbstractPolyVar}(v::Vector{PV}, x::AbstractMonomial)
    if v == x.vars
        # /!\ no copy done here for efficiency, do not mess up with vars
        w = v
        updatez = z -> z + x.z
    else
        w, maps = myunion([v, x.vars])
        updatez = z -> begin
            newz = zeros(Int, length(w))
            newz[maps[1]] += z
            newz[maps[2]] += x.z
            newz
        end
    end
    w, updatez
end

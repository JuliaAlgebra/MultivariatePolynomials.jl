function multiplyexistingvar{C}(v::Vector{PolyVar{C}}, x::PolyVar{C}, i::Int)
    updatez = z -> begin
        newz = copy(z)
        newz[i] += 1
        newz
    end
    # /!\ v not copied for efficiency, do not mess up with vars
    v, updatez
end
function insertvar{C}(v::Vector{PolyVar{C}}, x::PolyVar{C}, i::Int)
    n = length(v)
    I = 1:i-1
    J = i:n
    K = J+1
    w = Vector{PolyVar{C}}(n+1)
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

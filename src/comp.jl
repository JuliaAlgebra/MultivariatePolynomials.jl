import Base.==, Base.isless, Base.isapprox

# Comparison of variables

# TODO this should be in Base
function (==){T}(x::Vector{T}, y::Vector{T})
    if length(x) != length(y)
        false
    else
        #for (xi, yi) in zip(x, y)
        for i in 1:length(x)
            if x[i] != y[i]
                return false
            end
        end
        true
    end
end

# Technique: the higher catch the calls when it is rhs
# so p::PolyType == x::PolyVar -> x == p
(==)(p::PolyType, y) = y == p

# TODO equality should be between name ?
function (==){PV<:AbstractPolyVar}(x::PV, y::PV)
    x === y
end
(==)(α, x::AbstractPolyVar) = false
(==)(p::PolyType, x::AbstractPolyVar) = x == p

isless{PV<:AbstractPolyVar}(x::PV, y::PV) = isless(y.name, x.name)

# Comparison of monomials

# graded lex ordering
function mycomp{M<:AbstractMonomial}(x::M, y::M)
    degx = deg(x)
    degy = deg(y)
    if degx != degy
        degx - degy
    else
        i = j = 1
        # since they have the same degree,
        # if we get j > nvars(y), the rest in x.z should be zeros
        while i <= nvars(x) && j <= nvars(y)
            if x.vars[i] > y.vars[j]
                if x.z[i] == 0
                    i += 1
                else
                    return 1
                end
            elseif x.vars[i] < y.vars[j]
                if y.z[j] == 0
                    j += 1
                else
                    return -1
                end
            elseif x.z[i] != y.z[j]
                return x.z[i] - y.z[j]
            else
                i += 1
                j += 1
            end
        end
        0
    end
end

function (==){M<:AbstractMonomial}(x::M, y::M)
    mycomp(x, y) == 0
end
(==)(α, x::AbstractMonomial) = false
(==){M<:AbstractMonomial}(x::AbstractPolyVar, y::M) = M(x) == y
(==)(p::PolyType, x::AbstractMonomial) = x == p

# Comparison of MonomialVector
function (==){MV<:AbstractMonomialVector}(x::MV, y::MV)
    if length(x.Z) != length(y.Z)
        return false
    end
    allvars, maps = myunion([vars(x), vars(y)])
    # Should be sorted in the same order since the non-common
    # polyvar should have exponent 0
    for (a, b) in zip(x.Z, y.Z)
        A = zeros(length(allvars))
        B = zeros(length(allvars))
        A[maps[1]] = a
        B[maps[2]] = b
        if A != B
            return false
        end
    end
    return true
end
(==)(α, x::AbstractMonomialVector) = false
(==)(p::PolyType, x::AbstractMonomialVector) = false

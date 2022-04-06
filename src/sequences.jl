# Copied from TypedPolynomials
module Sequences

import MutableArithmetics
const MA = MutableArithmetics

export shortest_common_supersequence

"""
Given arrays (or strings) a and b, returns a table of size
(length(a) + 1) x (length(b) + 1)
Each entry i,j in the table is the length of the longest
common subsequence among a[1:(i-1)] and b[1:(j-1)]
Adapted from https://github.com/WestleyArgentum/Subsequences.jl
Copyright 2015 WestleyArentum and other contributors, released
under the MIT "Expat" License.
"""
function longest_common_subsequence_table(a, b)
    lengths = zeros(Int, length(a) + 1, length(b) + 1)
    for i in 1:length(a)
        ai = a[i]
        for j in 1:length(b)
            if ai == b[j]
                lengths[i+1, j+1] = lengths[i, j] + 1
            else
                lengths[i+1, j+1] = max(lengths[i+1, j], lengths[i, j+1])
            end
        end
    end
    lengths
end

"""
Given arrays (or strings) a and b, returns a table of size
(length(a) + 1) x (length(b) + 1)
Each entry i,j in the table is the length of the shortest
common supersequence among a[1:(i-1)] and b[1:(j-1)]
"""
function shortest_common_supersequence_table(a, b)
    table = longest_common_subsequence_table(a, b)
    for i in 1:size(table, 1)
        for j in 1:size(table, 2)
            table[i, j] = (i - 1) + (j - 1) - table[i, j]
        end
    end
    table
end

"""
Given arrays a and b, returns the shortest
common supersequence of a and b.
See: https://www.ics.uci.edu/~eppstein/161/960229.html
"""
function shortest_common_supersequence(a::Union{NTuple, AbstractVector}, b::Union{NTuple, AbstractVector})
    table = shortest_common_supersequence_table(a, b)
    result = promote_type(eltype(a), eltype(b))[]
    i = size(table, 1) - 1
    j = size(table, 2) - 1
    while i > 0 || j > 0
        if i == 0
            push!(result, b[j])
            j -= 1
        elseif j == 0
            push!(result, a[i])
            i -= 1
        elseif a[i] == b[j]
            push!(result, a[i])
            i -= 1
            j -= 1
        elseif table[i, j + 1] > table[i + 1, j]
            push!(result, b[j])
            j -= 1
        else
            push!(result, a[i])
            i -= 1
        end
    end
    reverse!(result)
    result
end

function shortest_common_supersequence(a::AbstractString, b::AbstractString)
    scs = shortest_common_supersequence(collect(a), collect(b))
    join(scs)
end

function mergesorted!(result, v1::AbstractArray, v2::AbstractArray, isless, combine, filter=x -> !iszero(x))
    i = 1
    i1 = 1
    i2 = 1
    while i1 <= length(v1) && i2 <= length(v2)
        x1 = v1[i1]
        x2 = v2[i2]
        if isless(x1, x2)
            if filter(x1)
                result[i] = x1
                i += 1
            end
            i1 += 1
        elseif isless(x2, x1)
            if filter(x2)
                result[i] = x2
                i += 1
            end
            i2 += 1
        else
            c = combine(x1, x2)
            if filter(c)
                result[i] = combine(x1, x2)
                i += 1
            end
            i1 += 1
            i2 += 1
        end
    end
    for j in i1:length(v1)
        if filter(v1[j])
            result[i] = v1[j]
            i += 1
        end
    end
    for j in i2:length(v2)
        if filter(v2[j])
            result[i] = v2[j]
            i += 1
        end
    end
    resize!(result, i - 1)
    result
end

function mergesorted(v1::AbstractArray, v2::AbstractArray, isless=Base.isless,
                     combine=Base.:+, filter=x -> !iszero(x))
    T = MA.promote_operation(combine, eltype(v1), eltype(v2))
    result = Vector{T}(undef, length(v1) + length(v2))
    mergesorted!(result, v1, v2, isless, combine, filter)
end

end

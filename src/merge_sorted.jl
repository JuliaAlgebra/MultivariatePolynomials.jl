import MergeSorted

function merge_sorted(a::AbstractVector, b::AbstractVector; kws...)
    return unique!(MergeSorted.mergesorted(a, b; kws...))
end

# TODO put it in MergeSorted.jl ?
merge_sorted(a::Tuple{}, b::Tuple{}; kws...) = tuple()
merge_sorted(a::Tuple, b::Tuple{}; kws...) = a
merge_sorted(a::Tuple{}, b::Tuple; kws...) = b
function merge_sorted(a::Tuple, b::Tuple; lt = isless, rev = false)
    x = first(a)
    y = first(b)
    return if x == y
        (x, merge_sorted(Base.tail(a), Base.tail(b); lt, rev)...)
    elseif xor(lt(x, y), rev)
        (x, merge_sorted(Base.tail(a), b; lt, rev)...)
    else
        (y, merge_sorted(a, Base.tail(b); lt, rev)...)
    end
end

function search_sorted_first(haystack::AbstractVector, needle; kws...)
    return searchsortedfirst(haystack, needle; kws...)
end
function search_sorted_first(haystack::Tuple, needle; kws...)
    return findfirst(isequal(needle), haystack)
end

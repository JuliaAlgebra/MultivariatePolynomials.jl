export Measure, zeta, ζ

type Measure{T}
  a::Vector{T}
  x::MonomialVector

  function Measure(a::Vector{T}, x::MonomialVector)
    if length(a) != length(x)
      error("There should be as many coefficient than monomials")
    end
    zeroidx = Int[]
    for (i,α) in enumerate(a)
      if iszero(α)
        push!(zeroidx, i)
      end
    end
    if !isempty(zeroidx)
      isnz = ones(Bool, length(a))
      isnz[zeroidx] = false
      nzidx = find(isnz)
      a = a[nzidx]
      x = x[nzidx]
    end
    new(a, x)
  end
end

Measure{T}(a::Vector{T}, x::MonomialVector) = Measure{T}(a, x)
function Measure(a::Vector, x::Vector)
    if length(a) != length(x)
        error("There should be as many coefficient than monomials")
    end
    perm, X = sortmonovec(x)
    Measure(a[perm], X)
end

function ζ{T}(v::Vector{T}, x::MonomialVector, varorder::Vector{PolyVar})
  Measure(T[m(v, varorder) for m in x], x)
end

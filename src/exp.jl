export Expectation

type Expectation{T}
  a::Vector{T}
  x::MonomialVector

  function Expectation(a::Vector{T}, x::MonomialVector)
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

Expectation{T}(a::Vector{T}, x::MonomialVector) = Expectation{T}(a, x)

function zeta{T}(v::Vector{T}, x::MonomialVector, varorder::Vector{PolyVar})
  Expectation(T[m(v, varorder) for m in x], x)
end

ζ(v::Vector, x::MonomialVector, varorder::Vector{PolyVar}) = zeta(v, x, varorder)

function dot(pexp::Expectation, p::TermContainer)
  i = 1
  s = 0
  for t in p
    while t.x != pexp.x[i] && i <= length(pexp.x)
      i += 1
    end
    if i > length(pexp.x)
      error("The polynomial $p has a monomial for which the expectation is not known in $pexp")
    end
    s += pexp.a * t.α
    i += 1
  end
  s
end
dot(p::TermContainer, pexp::Expectation) = dot(pexp, p)
dot(pexp::Expectation, p::PolyType) = dot(pexp, TermContainer(p))
dot(p::PolyType, pexp::Expectation) = dot(pexp, TermContainer(p))
(pexp::Expectation)(p::PolyType) = dot(pexp, p)

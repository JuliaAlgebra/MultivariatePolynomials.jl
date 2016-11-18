function dot(m::Measure, p::TermContainer)
  i = 1
  s = 0
  for t in p
    while i <= length(m.x) && t.x != m.x[i]
      i += 1
    end
    if i > length(m.x)
      error("The polynomial $p has a monomial for which the expectation is not known in $m")
    end
    s += m.a[i] * t.α
    i += 1
  end
  s
end
dot(p::TermContainer, m::Measure) = dot(m, p)
dot(m::Measure, p::PolyType) = dot(m, TermContainer(p))
dot(p::PolyType, m::Measure) = dot(m, TermContainer(p))
exp(m::Measure, p::PolyType) = dot(m, p)
exp(p::PolyType, m::Measure) = dot(m, p)

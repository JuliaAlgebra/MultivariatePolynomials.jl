export expectation

function _dot{C}(m::Measure{C}, p::TermContainer{C}, f)
    i = 1
    s = 0
    for t in p
        while i <= length(m.x) && t.x != m.x[i]
            i += 1
        end
        if i > length(m.x)
            error("The polynomial $p has a nonzero term $t with monomial $(t.x) for which the expectation is not known in $m")
        end
        s += f(m.a[i], t.α)
        i += 1
    end
    s
end
dot(m::Measure, p::TermContainer) = _dot(m, p, (*))
dot(p::TermContainer, m::Measure) = _dot(m, p, (a, b) -> b * a)

dot(m::Measure, p::PolyType) = dot(m, TermContainer(p))
dot(p::PolyType, m::Measure) = dot(TermContainer(p), m)

expectation(m::Measure, p::PolyType) = dot(m, p)
expectation(p::PolyType, m::Measure) = dot(p, m)

export expectation

function _dot(m::AbstractMeasure, p::AbstractTermContainer, f)
    i = 1
    s = 0
    for t in p
        while i <= length(m.x) && t.x != m.x[i]
            i += 1
        end
        if i > length(m.x)
            error("The polynomial $p has a monomial for which the expectation is not known in $m")
        end
        s += f(m.a[i], t.α)
        i += 1
    end
    s
end
dot(m::AbstractMeasure, p::AbstractTermContainer) = _dot(m, p, (*))
dot(p::AbstractTermContainer, m::AbstractMeasure) = _dot(m, p, (a, b) -> b * a)

tcfor(::Measure) = TermContainer
tcfor(::NCMeasure) = NCTermContainer
dot(m::AbstractMeasure, p::PolyType) = dot(m, tcfor(m)(p))
dot(p::PolyType, m::AbstractMeasure) = dot(tcfor(m)(p), m)

expectation(m::AbstractMeasure, p::PolyType) = dot(m, p)
expectation(p::PolyType, m::AbstractMeasure) = dot(p, m)

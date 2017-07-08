export expectation

function _dot(m::Measure, p::APL, f)
    i = 1
    s = 0
    for t in terms(p)
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
Base.dot(m::Measure, p::APL) = _dot(m, p, (*))
Base.dot(p::APL, m::Measure) = _dot(m, p, (a, b) -> b * a)

expectation(m::Measure, p::APL) = dot(m, p)
expectation(p::APL, m::Measure) = dot(p, m)

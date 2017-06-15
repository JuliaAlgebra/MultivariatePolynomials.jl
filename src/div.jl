# _div(a, b) assumes that b divides a
function _div(m1::Monomial{true}, m2::Monomial{true})
    w, updatez = multdivmono(m1.vars, m2, -)
    Monomial{true}(w, updatez(m1.z))
end
function _div(t::Term, m::Monomial)
    t1.α * _div(monomial(t), m)
end
function _div(t1::Term, t2::Term)
    (t1.α / t2.α) * _div(monomial(t1), monomial(t2))
end

proddiff(x, y) = x*y - x*y

function Base.divrem{C, T, S}(f::Polynomial{C, T}, g::Polynomial{C, S})
    rf = Polynomial{C, promote_op(proddiff, T, S)}(f)
    q = r = zero(f - g / 2)
    lt = leadingterm(g)
    rg = removeleadingterm(g)
    lm = monomial(lt)
    while !iszero(rf)
        ltf = leadingterm(rf)
        if divides(lm, ltf)
            qt = _div(ltf, lt)
            q += qt
            rf = removeleadingterm(rf) - qt * rg
        elseif lm > monomial(ltf)
            # Since the monomials are sorted in decreasing order,
            # lm is larger than all of them hence it cannot divide any of them
            r += rf
            break
        else
            r += ltf
            rf = removeleadingterm(rf)
        end
    end
    q, r
end
Base.div(f::PolyType, g::PolyType) = divrem(f, g)[1]
Base.rem(f::PolyType, g::PolyType) = divrem(f, g)[2]

function evalmap(vars, x::Vector, varorder::Vector{PolyVar})
    vals = Any[var for var in vars]
    for (i, var) in enumerate(varorder)
        j = findfirst(vars, var)
        # If i == 0, that means that the variable is not present
        # so it is ignored
        if j > 0
            vals[j] = x[i]
        end
    end
    vals
end

function termeval(t::Term, vals::Vector)
    val = t.α
    for i in 1:length(vals)
        if t.x.z[i] > 0
            val *= vals[i]^t.x.z[i]
        end
    end
    val
end

function (m::PolyVar)(x::Vector, varorder)
    Term(m)(x, varorder)
end

function (m::Monomial)(x::Vector, varorder)
    Term(m)(x, varorder)
end

function (t::Term)(x::Vector, varorder)
    vals = evalmap(vars(t), x, varorder)
    termeval(t, vals)
end

function (p::VecPolynomial)(x::Vector, varorder)
    vals = evalmap(vars(p), x, varorder)
    q = zero(p)
    for t in p
        q += termeval(t, vals)
    end
    q
end

function (p::MatPolynomial)(x::Vector, varorder)
    VecPolynomial(p)(x, varorder)
end

function (q::RationalPoly)(x::Vector, varorder)
    q.num(x, varorder) / q.den(x, varorder)
end

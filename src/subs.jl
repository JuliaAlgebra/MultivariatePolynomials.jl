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

function monoeval(z::Vector{Int}, vals::Vector)
    @assert length(z) == length(vals)
    @assert !isempty(z)
    val = vals[1]^z[1]
    for i in 2:length(vals)
        if z[i] > 0
            val *= vals[i]^z[i]
        end
    end
    val
end

function (m::PolyVar)(x::Vector, varorder)
    Monomial(m)(x, varorder)
end

function (m::Monomial)(x::Vector, varorder)
    vals = evalmap(vars(m), x, varorder)
    monoeval(m.z, vals)
end

function (t::Term)(x::Vector, varorder)
    vals = evalmap(vars(t), x, varorder)
    t.α * monoeval(t.x.z, vals)
end

function (p::VecPolynomial)(x::Vector, varorder)
    vals = evalmap(vars(p), x, varorder)
    q = zero(p)
    for i in 1:length(p)
        q += p.a[i] * monoeval(p.x.Z[i], vals)
    end
    q
end

function (p::MatPolynomial)(x::Vector, varorder)
    VecPolynomial(p)(x, varorder)
end

function (q::RationalPoly)(x::Vector, varorder)
    q.num(x, varorder) / q.den(x, varorder)
end

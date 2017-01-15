# eval replace all variables by a new value and subs replace some of them.
# we use different function for that because of type inferability.
# with eval, Julia knows that all PolyVar will be replaced by values so it can do better inference.
export subs

function fillmap!(vals, vars, x, varorder)
    for (i, var) in enumerate(varorder)
        j = findfirst(vars, var)
        # If i == 0, that means that the variable is not present
        # so it is ignored
        if j > 0
          vals[j] = x[i]
        end
    end
end

function evalmap{T}(vars, x::Vector{T}, varorder::Vector{PolyVar})
    if vars == varorder
        x
    else
        # Every variable will be replaced by some value of type T
        vals = Vector{T}(length(vars))
        fillmap!(vals, vars, x, varorder)
        for i in 1:length(vals)
          @assert isdefined(vals, i) "Variable $(vars[i]) was not assigned a value"
        end
        vals
    end
end

function subsmap{T}(vars, x::Vector{T}, varorder::Vector{PolyVar})
    if vars == varorder
        x
    else
        # Some variable may not be replaced
        vals = Vector{promote_type(T, PolyVar)}(length(vars))
        copy!(vals, vars)
        fillmap!(vals, vars, x, varorder)
        vals
    end
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
    sum(i -> p.a[i] * monoeval(p.x.Z[i], vals), 1:length(p))
end

function subs(p::VecPolynomial, x::Vector, varorder)
    vals = subsmap(vars(p), x, varorder)
    sum(i -> p.a[i] * monoeval(p.x.Z[i], vals), 1:length(p))
end

function (p::MatPolynomial)(x::Vector, varorder)
    VecPolynomial(p)(x, varorder)
end

function (q::RationalPoly)(x::Vector, varorder)
    q.num(x, varorder) / q.den(x, varorder)
end

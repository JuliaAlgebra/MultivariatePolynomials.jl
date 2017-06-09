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

function evalmap{C, T}(vars, x::Vector{T}, varorder::Vector{PolyVar{C}})
    if vars == varorder
        x
    else
        # Every variable will be replaced by some value of type T
        vals = Vector{T}(length(vars))
        fillmap!(vals, vars, x, varorder)
        for i in 1:length(vals)
            @assert isassigned(vals, i) "Variable $(vars[i]) was not assigned a value"
        end
        vals
    end
end

function subsmap{C, T}(vars, x::Vector{T}, varorder::Vector{PolyVar{C}})
    if vars == varorder
        x
    else
        # Some variable may not be replaced
        vals = Vector{promote_type(T, PolyVar{C})}(length(vars))
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

function (p::Polynomial{C, T}){C, T, S}(x::Vector{S}, varorder)
    vals = evalmap(vars(p), x, varorder)
    # I need to check for izero otherwise I get : ArgumentError: reducing over an empty collection is not allowed
    iszero(p) ? zero(Base.promote_op(*, S, T)) : sum(i -> p.a[i] * monoeval(p.x.Z[i], vals), 1:length(p))
end

function subs{C, T, S}(p::Polynomial{C, T}, x::Vector{S}, varorder)
    Tin = S <: PolyType ? S : eltype(S)
    Tout = Base.promote_op(*, T, Tin)
    vals = subsmap(vars(p), x, varorder)
    # I need to check for izero otherwise I get : ArgumentError: reducing over an empty collection is not allowed
    if iszero(p)
        zero(Polynomial{C, Tout})
    else
        Polynomial{C, Tout}(sum(i -> p.a[i] * monoeval(p.x.Z[i], vals), 1:length(p)))
    end
end

(p::MatPolynomial)(x::Vector, varorder) = Polynomial(p)(x, varorder)
subs(p::MatPolynomial, x::Vector, varorder) = subs(Polynomial(p), x, varorder)

function (q::RationalPoly)(x::Vector, varorder)
    q.num(x, varorder) / q.den(x, varorder)
end
subs(q::RationalPoly, x::Vector, varorder) = subs(q.num, x, varorder) / subs(q.den, x, varorder)

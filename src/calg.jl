# Multiple Dispatch strategy:
# Polytype op Any
# for each concrete type T
# Any op T # ambig with Polytype op Any
# PolyType op T # breaks ambig

import Base.(^), Base.dot

import Base.transpose
Base.transpose(p::PolyType) = p

# In Base/intfuncs.jl, x^p returns zero(x) when p == 0
# Since one(PolyVar) and one(Monomial) do not return
# a PolyVar and a Monomial, this results in type instability
# Defining the specific methods solve this problem and also make
# them a lot faster
^{C}(x::PolyVar{C}, i::Int) = Monomial{C}([x], [i])
^(x::Monomial{true}, i::Int) = Monomial{true}(x.vars, i*x.z)

# Product between PolyVar and Monomial -> Monomial
function (*)(x::PolyVar{true}, y::PolyVar{true})
    if x === y
        Monomial{true}([x], [2])
    else
        Monomial{true}(x > y ? [x,y] : [y,x], [1,1])
    end
end
function multiplyvar(v::Vector{PolyVar{true}}, x::PolyVar{true})
    i = findfirst(w->w <= x, v)
    if i > 0 && v[i] == x
        multiplyexistingvar(v, x, i)
    else
        insertvar(v, x, i == 0 ? length(v)+1 : i)
    end
end
function (*)(x::PolyVar{true}, y::Monomial{true})
    w, updatez = multiplyvar(y.vars, x)
    Monomial{true}(w, updatez(y.z))
end
function (*)(x::PolyVar{true}, y::MonomialVector{true})
    w, updatez = multiplyvar(y.vars, x)
    MonomialVector{true}(w, updatez.(y.Z))
end
function multdivmono(v::Vector{PolyVar{true}}, x::Monomial{true}, op)
    if v == x.vars
        # /!\ no copy done here for efficiency, do not mess up with vars
        w = v
        updatez = z -> op(z, x.z)
    else
        w, maps = myunion([v, x.vars])
        updatez = z -> begin
            newz = zeros(Int, length(w))
            newz[maps[1]] = op(newz[maps[1]], z)
            newz[maps[2]] = op(newz[maps[2]], x.z)
            newz
        end
    end
    w, updatez
end
function (*)(x::Monomial{true}, y::Monomial{true})
    w, updatez = multdivmono(y.vars, x, +)
    Monomial{true}(w, updatez(y.z))
end
function (*)(x::Monomial{true}, y::MonomialVector{true})
    w, updatez = multdivmono(y.vars, x, +)
    MonomialVector{true}(w, updatez.(y.Z))
end
(*)(x::Monomial{true}, y::PolyVar{true}) = y * x

# non-PolyType * PolyType: specific methods for speed
*(p::PolyType, α) = α * p

*{C}(α, x::PolyVar{C}) = Term(α, Monomial{C}(x))
*(α, x::Monomial)      = Term(α, x)
*(α, p::MatPolynomial) = α * Polynomial(p)
*(α, x::Term)    = Term(α*x.α, x.x)
*(α, p::Polynomial) = Polynomial(α*p.a, p.x)

# Reverse order to avoid abiguïty with above 5 specific methods
*(p::PolyType, x::PolyVar) = x * p
*(p::PolyType, x::Monomial) = x * p
*(p::PolyType, x::MatPolynomial) = x * Polynomial(p)
# The three above are mapped to one of the two below
*(p::PolyType, q::Term) = TermContainer(p) * q
*(p::PolyType, q::Polynomial) = TermContainer(p) * q

# I do not want to cast x to TermContainer because that would force the promotion of eltype(q) with Int
function *{S<:Union{PolyVar,Monomial}}(x::S, t::Term)
    Term(t.α, x*t.x)
end
function *{S<:Union{PolyVar,Monomial},T}(x::S, p::Polynomial{T})
    # /!\ No copy of a is done
    Polynomial{T}(p.a, x*p.x)
end

# TermContainer * TermContainer
*(x::Term, y::Term) = Term(x.α*y.α, x.x*y.x)
*(p::Polynomial, t::Term) = t * p
function *(t::Term, p::Polynomial)
    if iszero(t)
        zero(p)
    else
        n = length(p)
        allvars, maps = myunion([t.x.vars, p.x.vars])
        nvars = length(allvars)
        Z = [zeros(Int, nvars) for i in 1:n]
        for i in 1:n
            Z[i][maps[1]] = t.x.z
            Z[i][maps[2]] += p.x.Z[i]
        end
        Polynomial(t.α * p.a, MonomialVector(allvars, Z))
    end
end
function *(p::Polynomial, q::Polynomial)
    if iszero(p)
        zero(q)
    elseif iszero(q)
        zero(p)
    else
        samevars = vars(p) == vars(q)
        if samevars
            allvars = vars(p)
        else
            allvars, maps = myunion([vars(p), vars(q)])
        end
        N = length(p)*length(q)
        Z = Vector{Vector{Int}}(N)
        T = typeof(p.a[1]*q.a[1])
        a = Vector{T}(N)
        i = 0
        for u in p
            for v in q
                if samevars
                    z = u.x.z + v.x.z
                else
                    z = zeros(Int, length(allvars))
                    z[maps[1]] += u.x.z
                    z[maps[2]] += v.x.z
                end
                i += 1
                Z[i] = z
                a[i] = u.α * v.α
            end
        end
        vecpolynomialclean(allvars, a, Z)
    end
end

dot(p::PolyType, q::PolyType) = p * q
dot(α, p::PolyType) = α * p
dot(p::PolyType, α) = p * α

myminivect{T}(x::T, y::T) = [x, y]
function myminivect{S,T}(x::S, y::T)
    U = promote_type(S, T)
    [U(x), U(y)]
end

function (+){C}(x::Term{C}, y::Term{C})
    if x.x == y.x
        Polynomial{C}([x.α+y.α], [x.x])
    elseif x.x > y.x
        Polynomial{C}(myminivect(x.α,y.α), [x.x,y.x])
    else
        Polynomial{C}(myminivect(y.α,x.α), [y.x,x.x])
    end
end

function (-){C}(x::Term{C}, y::Term{C})
    if x.x == y.x
        Polynomial{C}([x.α-y.α], [x.x])
    elseif x.x > y.x
        Polynomial{C}(myminivect(x.α,-y.α), [x.x,y.x])
    else
        Polynomial{C}(myminivect(-y.α,x.α), [y.x,x.x])
    end
end

(+){S<:Union{PolyVar,Monomial},T<:Union{PolyVar,Monomial}}(x::S, y::T) = Term(x) + Term(y)
(-){S<:Union{PolyVar,Monomial},T<:Union{PolyVar,Monomial}}(x::S, y::T) = Term(x) - Term(y)

function plusorminus{C, S, T}(p::TermContainer{C, S}, q::TermContainer{C, T}, isplus)
    varsvec = [vars(p), vars(q)]
    allvars, maps = myunion(varsvec)
    nvars = length(allvars)
    U = promote_type(S, T)
    a = Vector{U}()
    Z = Vector{Vector{Int}}()
    i = j = 1
    while i <= length(p) || j <= length(q)
        z = zeros(Int, nvars)
        if j > length(q) || (i <= length(p) && p[i].x > q[j].x)
            t = p[i]
            z[maps[1]] = t.x.z
            α = t.α
            i += 1
        elseif i > length(p) || q[j].x > p[i].x
            t = q[j]
            z[maps[2]] = t.x.z
            α = isplus ? t.α : -t.α
            j += 1
        else
            t = p[i]
            z[maps[1]] = t.x.z
            α = t.α
            s = q[j]
            α += isplus ? s.α : -s.α
            i += 1
            j += 1
        end
        push!(a, α)
        push!(Z, z)
    end

    Polynomial(a, MonomialVector{C}(allvars, Z))
end


(+){C}(x::TermContainer{C}, y::TermContainer{C}) = plusorminus(x, y, true)
(-){C}(x::TermContainer{C}, y::TermContainer{C}) = plusorminus(x, y, false)
(+){C, T, S<:Union{Monomial,PolyVar}}(x::TermContainer{C, T}, y::S) = x + Term{C, T}(y)
(+){C, T, S<:Union{Monomial,PolyVar}}(x::S, y::TermContainer{C, T}) = Term{C, T}(x) + y

(+)(x::TermContainer, y::MatPolynomial) = x + Polynomial(y)
(+)(x::MatPolynomial, y::TermContainer) = Polynomial(x) + y
(+)(x::MatPolynomial, y::MatPolynomial) = Polynomial(x) + Polynomial(y)
(-)(x::TermContainer, y::MatPolynomial) = x - Polynomial(y)
(-)(x::MatPolynomial, y::TermContainer) = Polynomial(x) - y
(-)(x::MatPolynomial, y::MatPolynomial) = Polynomial(x) - Polynomial(y)

(-){S<:Union{Monomial,PolyVar},T}(x::TermContainer{T}, y::S) = x - Term{T}(y)
(-){S<:Union{Monomial,PolyVar},T}(x::S, y::TermContainer{T}) = Term{T}(x) - y

# Avoid adding a zero constant that might artificially increase the Newton polytope
# Need to add Polynomial conversion for type stability
(+){C}(x::PolyType{C}, y) = iszero(y) ? Polynomial(x) : x + Term{C}(y)
(+){C}(x, y::PolyType{C}) = iszero(x) ? Polynomial(y) : Term{C}(x) + y
(-){C}(x::PolyType{C}, y) = iszero(y) ? Polynomial(x) : x - Term{C}(y)
(-){C}(x, y::PolyType{C}) = iszero(x) ? Polynomial(-y) : Term{C}(x) - y

(-){C}(x::PolyVar{C}) = Term(-1, Monomial{C}(x))
(-)(x::Monomial) = Term(-1, x)
(-)(t::Term) = Term(-t.α, t.x)
(-)(p::Polynomial) = Polynomial(-p.a, p.x)

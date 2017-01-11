# Multiple Dispatch strategy:
# Polytype op Any
# for each concrete type T
# Any op T # ambig with Polytype op Any
# PolyType op T # breaks ambig

import Base.(^), Base.dot, Base.(.+), Base.(.-), Base.(.*), Base.(./), Base.(.^)

(.+)(p::PolyType, α) = p+α
(.+)(α, p::PolyType) = α+p
(.+)(p::PolyType, q::PolyType) = p+q
(.-)(p::PolyType, α) = p-α
(.-)(α, p::PolyType) = α-p
(.-)(p::PolyType, q::PolyType) = p-q
(.*)(p::PolyType, α) = p*α
(.*)(α, p::PolyType) = α*p
(.*)(p::PolyType, q::PolyType) = p*q
(./)(p::PolyType, α) = p/α
(./)(α, p::PolyType) = α/p
(./)(p::PolyType, q::PolyType) = p/q
(.^)(p::PolyType, i::Int) = p^i

import Base.transpose
Base.transpose(p::PolyType) = p

# In Base/intfuncs.jl, x^p returns zero(x) when p == 0
# Since one(PolyVar) and one(Monomial) do not return
# a PolyVar and a Monomial, this results in type instability
# Defining the specific methods solve this problem and also make
# them a lot faster
^(x::PolyVar, i::Int) = Monomial([x], [i])
^(x::Monomial, i::Int) = Monomial(x.vars, i*x.z)

# Product between PolyVar and Monomial -> Monomial
function (*)(x::PolyVar, y::PolyVar)
    if x === y
        Monomial([x], [2])
    else
        Monomial(x > y ? [x,y] : [y,x], [1,1])
    end
end
function multiplyvar(v::Vector{PolyVar}, x::PolyVar)
    i = findfirst(w->w <= x, v)
    if i > 0 && v[i] == x
        multiplyexistingvar(v, x, i)
    else
        insertvar(v, x, i == 0 ? length(v)+1 : i)
    end
end
function (*)(x::PolyVar, y::Monomial)
    w, updatez = multiplyvar(y.vars, x)
    Monomial(w, updatez(y.z))
end
function (*)(x::PolyVar, y::MonomialVector)
    w, updatez = multiplyvar(y.vars, x)
    MonomialVector(w, updatez.(y.Z))
end
function multiplymono(v::Vector{PolyVar}, x::Monomial)
    if v == x.vars
        # /!\ no copy done here for efficiency, do not mess up with vars
        w = v
        updatez = z -> z + x.z
    else
        w, maps = myunion([v, x.vars])
        updatez = z -> begin
            newz = zeros(Int, length(w))
            newz[maps[1]] += z
            newz[maps[2]] += x.z
            newz
        end
    end
    w, updatez
end
function (*)(x::Monomial, y::Monomial)
    w, updatez = multiplymono(y.vars, x)
    Monomial(w, updatez(y.z))
end
function (*)(x::Monomial, y::MonomialVector)
    w, updatez = multiplymono(y.vars, x)
    MonomialVector(w, updatez.(y.Z))
end
(*)(x::Monomial, y::PolyVar) = y * x

# non-PolyType * PolyType: specific methods for speed
*(p::PolyType, α) = α * p

*(α, x::PolyVar)       = Term(α, Monomial(x))
*(α, x::Monomial)      = Term(α, x)
*(α, p::MatPolynomial) = α * VecPolynomial(p)
*{T}(α, x::Term{T})    = Term(T(α)*x.α, x.x)
*(α, p::VecPolynomial) = VecPolynomial(α*p.a, p.x)

# Reverse order to avoid abiguïty with above 5 specific methods
*(p::PolyType, x::PolyVar) = x * p
*(p::PolyType, x::Monomial) = x * p
*(p::PolyType, x::MatPolynomial) = x * VecPolynomial(p)
# The three above are mapped to one of the two below
*(p::PolyType, q::Term) = TermContainer(p) * q
*(p::PolyType, q::VecPolynomial) = TermContainer(p) * q

# I do not want to cast x to TermContainer because that would force the promotion of eltype(q) with Int
function *{S<:Union{PolyVar,Monomial}}(x::S, t::Term)
    Term(t.α, x*t.x)
end
function *{S<:Union{PolyVar,Monomial},T}(x::S, p::VecPolynomial{T})
    # /!\ No copy of a is done
    VecPolynomial{T}(p.a, x*p.x)
end

# TermContainer * TermContainer
*(x::Term, y::Term) = Term(x.α*y.α, x.x*y.x)
*(p::VecPolynomial, t::Term) = t * p
function *(t::Term, p::VecPolynomial)
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
        VecPolynomial(t.α * p.a, MonomialVector(allvars, Z))
    end
end
function *(p::VecPolynomial, q::VecPolynomial)
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

function (+)(x::Term, y::Term)
    if x.x == y.x
        VecPolynomial([x.α+y.α], [x.x])
    elseif x.x > y.x
        VecPolynomial(myminivect(x.α,y.α), [x.x,y.x])
    else
        VecPolynomial(myminivect(y.α,x.α), [y.x,x.x])
    end
end

function (-)(x::Term, y::Term)
    if x.x == y.x
        VecPolynomial([x.α-y.α], [x.x])
    elseif x.x > y.x
        VecPolynomial(myminivect(x.α,-y.α), [x.x,y.x])
    else
        VecPolynomial(myminivect(-y.α,x.α), [y.x,x.x])
    end
end

(+){S<:Union{PolyVar,Monomial},T<:Union{PolyVar,Monomial}}(x::S, y::T) = Term(x) + Term(y)
(-){S<:Union{PolyVar,Monomial},T<:Union{PolyVar,Monomial}}(x::S, y::T) = Term(x) - Term(y)

function plusorminus{S,T}(p::TermContainer{S}, q::TermContainer{T}, isplus)
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

    VecPolynomial(a, MonomialVector(allvars, Z))
end


(+)(x::TermContainer, y::TermContainer) = plusorminus(x, y, true)
(-)(x::TermContainer, y::TermContainer) = plusorminus(x, y, false)
(+){S<:Union{Monomial,PolyVar},T}(x::TermContainer{T}, y::S) = x + Term{T}(y)
(+){S<:Union{Monomial,PolyVar},T}(x::S, y::TermContainer{T}) = Term{T}(x) + y

(+)(x::TermContainer, y::MatPolynomial) = x + VecPolynomial(y)
(+)(x::MatPolynomial, y::TermContainer) = VecPolynomial(x) + y
(+)(x::MatPolynomial, y::MatPolynomial) = VecPolynomial(x) + VecPolynomial(y)
(-)(x::TermContainer, y::MatPolynomial) = x - VecPolynomial(y)
(-)(x::MatPolynomial, y::TermContainer) = VecPolynomial(x) - y
(-)(x::MatPolynomial, y::MatPolynomial) = VecPolynomial(x) - VecPolynomial(y)

iszero{T}(x::T) = x == zero(T)
iszero(t::Term) = iszero(t.α)
iszero(p::VecPolynomial) = isempty(p.x)
iszero(p::MatPolynomial) = isempty(VecPolynomial(p).x)

(-){S<:Union{Monomial,PolyVar},T}(x::TermContainer{T}, y::S) = x - Term{T}(y)
(-){S<:Union{Monomial,PolyVar},T}(x::S, y::TermContainer{T}) = Term{T}(x) - y

# Avoid adding a zero constant that might artificially increase the Newton polytope
# Need to add VecPolynomial conversion for type stability
(+)(x::PolyType, y) = iszero(y) ? VecPolynomial(x) : x + Term(y)
(+)(x, y::PolyType) = iszero(x) ? VecPolynomial(y) : Term(x) + y
(-)(x::PolyType, y) = iszero(y) ? VecPolynomial(x) : x - Term(y)
(-)(x, y::PolyType) = iszero(x) ? VecPolynomial(-y) : Term(x) - y

(-)(x::PolyVar) = Term(-1, Monomial(x))
(-)(x::Monomial) = Term(-1, x)
(-)(t::Term) = Term(-t.α, t.x)
(-)(p::VecPolynomial) = VecPolynomial(-p.a, p.x)

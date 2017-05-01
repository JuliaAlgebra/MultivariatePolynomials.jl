import Base.==, Base.isless, Base.isapprox
export isapproxzero

# iszero is only available in Julia v0.6
if isdefined(Base, :iszero)
    import Base.iszero
else
    iszero{T}(x::T) = x == zero(T)
end
iszero(t::Term) = iszero(t.α)
iszero(p::Polynomial) = isempty(p)
iszero(p::MatPolynomial) = isempty(Polynomial(p))

# TODO This should be in Base with T instead of PolyVar{C}.
# See https://github.com/blegat/MultivariatePolynomials.jl/issues/3
function (==){C}(x::Vector{PolyVar{C}}, y::Vector{PolyVar{C}})
    if length(x) != length(y)
        false
    else
        #for (xi, yi) in zip(x, y)
        for i in 1:length(x)
            if x[i] != y[i]
                return false
            end
        end
        true
    end
end

# Technique: the higher catch the calls when it is rhs
# so p::PolyType == x::PolyVar -> x == p
(==)(p::PolyType, y) = y == p

# Comparison of PolyVar

function (==){C}(x::PolyVar{C}, y::PolyVar{C})
    x.id == y.id
end
(==)(α, x::PolyVar) = false
(==){C}(p::PolyType{C}, x::PolyVar{C}) = x == p

isless{C}(x::PolyVar{C}, y::PolyVar{C}) = isless(y.id, x.id)

# Comparison of Monomial

# graded lex ordering
function mycomp{C}(x::Monomial{C}, y::Monomial{C})
    degx = deg(x)
    degy = deg(y)
    if degx != degy
        degx - degy
    else
        i = j = 1
        # since they have the same degree,
        # if we get j > nvars(y), the rest in x.z should be zeros
        while i <= nvars(x) && j <= nvars(y)
            if x.vars[i] > y.vars[j]
                if x.z[i] == 0
                    i += 1
                else
                    return 1
                end
            elseif x.vars[i] < y.vars[j]
                if y.z[j] == 0
                    j += 1
                else
                    return -1
                end
            elseif x.z[i] != y.z[j]
                return x.z[i] - y.z[j]
            else
                i += 1
                j += 1
            end
        end
        0
    end
end

function (==){C}(x::Monomial{C}, y::Monomial{C})
    mycomp(x, y) == 0
end
(==)(α, x::Monomial) = isconstant(x) && α == 1
(==){C}(x::PolyVar{C}, y::Monomial{C}) = Monomial{C}(x) == y
(==){C}(p::PolyType{C}, x::Monomial{C}) = x == p

# graded lex ordering
function isless{C}(x::Monomial{C}, y::Monomial{C})
    mycomp(x, y) < 0
end
isless{C}(x::Monomial{C}, y::PolyVar{C}) = isless(x, Monomial{C}(y))
isless{C}(x::PolyVar{C}, y::Monomial{C}) = isless(Monomial{C}(x), y)

# Comparison of MonomialVector
function (==){C}(x::MonomialVector{C}, y::MonomialVector{C})
    if length(x.Z) != length(y.Z)
        return false
    end
    allvars, maps = myunion([vars(x), vars(y)])
    # Should be sorted in the same order since the non-common
    # polyvar should have exponent 0
    for (a, b) in zip(x.Z, y.Z)
        A = zeros(length(allvars))
        B = zeros(length(allvars))
        A[maps[1]] = a
        B[maps[2]] = b
        if A != B
            return false
        end
    end
    return true
end
(==)(α, x::MonomialVector) = false
(==){C}(p::PolyType{C}, x::MonomialVector{C}) = false

# Comparison of Term
function isapproxzero(x; ztol::Real=1e-6)
    -ztol < x < ztol
end

# See https://github.com/blegat/MultivariatePolynomials.jl/issues/22
(==)(α::Void, x::TermType) = false
(==)(α::Dict, x::TermType) = false
(==)(x::TermType, α::Dict) = false
(==){C}(y, p::TermType{C}) = TermContainer{C}(y) == p
(==)(y::PolyType, p::TermContainer) = TermContainer(y) == p

function (==){C}(s::Term{C}, t::Term{C})
    (s.α == t.α) && (iszero(s.α) || s.x == t.x)
end
function (==){C}(t::Term{C}, p::Polynomial{C})
    if isempty(p.a)
        iszero(t.α)
    else
        length(p.a) == 1 && p.a[1] == t.α && p.x[1] == t.x
    end
end
(==){C}(p::Polynomial{C}, t::Term{C}) = t == p
function (==){C}(p::Polynomial{C}, q::Polynomial{C})
    # terms should be sorted and without zeros
    if length(p) != length(q)
        return false
    end
    for i in 1:length(p)
        if p.x[i] != q.x[i]
            # There should not be zero terms
            @assert p.a[i] != 0
            @assert q.a[i] != 0
            return false
        end
        if p.a[i] != q.a[i]
            return false
        end
    end
    return true
end

(==)(p::RationalPoly, q::RationalPoly) = p.num*q.den == q.num*p.den
# Solve ambiguity with (::PolyType, ::Any)
(==)(p::PolyType, q::RationalPoly) = p*q.den == q.num
(==)(p, q::RationalPoly) = p*q.den == q.num
# IJulia output, see https://github.com/blegat/MultivariatePolynomials.jl/issues/22
(==)(α::Void, x::RationalPoly) = false
(==)(α::Dict, x::RationalPoly) = false

(==)(p::TermContainer, q::MatPolynomial) = p == TermContainer(q)
(==)(p::MatPolynomial, q::MatPolynomial) = iszero(p - q)

function grlex(x::Vector{Int}, y::Vector{Int})
    @assert length(x) == length(y)
    degx = sum(x)
    degy = sum(y)
    if degx != degy
        degx < degy
    else
        for (a, b) in zip(x, y)
            if a < b
                return true
            elseif a > b
                return false
            end
        end
        false
    end
end

function isapproxzero(p::Polynomial; ztol::Real=1e-6)
    isapprox(p, zero(p), ztol=ztol)
end

function isapproxzero(p::RationalPoly; ztol::Real=1e-6)
    isapproxzero(p.num, ztol=ztol)
end

function isapprox{C, S, T}(p::Polynomial{C, S}, q::Polynomial{C, T}; rtol::Real=Base.rtoldefault(S, T), atol::Real=0, ztol::Real=1e-6)
    i = j = 1
    while i <= length(p.x) || j <= length(q.x)
        lhs, rhs = 0, 0
        if i > length(p.x) || (j <= length(q.x) && q.x[j] > p.x[i])
            if !isapproxzero(q.a[j], ztol=ztol)
                return false
            end
            j += 1
        elseif j > length(q.x) || p.x[i] > q.x[j]
            if !isapproxzero(p.a[i], ztol=ztol)
                return false
            end
            i += 1
        else
            if !isapprox(p.a[i], q.a[j], rtol=rtol, atol=atol)
                return false
            end
            i += 1
            j += 1
        end
    end
    true
end

function isapprox{C, S, T}(s::Term{C, S}, t::Term{C, T}; rtol::Real=Base.rtoldefault(S, T), atol::Real=0, ztol::Real=1e-6)
    s.x == t.x && isapprox(s.α, t.α, rtol=rtol, atol=atol)
end

function isapprox(p::MatPolynomial, q::MatPolynomial)
    p.x == q.x && isapprox(p.Q, q.Q)
end

function permcomp(f, m)
    picked = IntSet()
    for i in 1:m
        k = 0
        for j in 1:m
            if !(j in picked) && f(i, j)
                k = j
                break
            end
        end
        if k == 0
            return false
        end
        push!(picked, k)
    end
    true
end

function isapprox{C, S, T}(p::SOSDecomposition{C, S}, q::SOSDecomposition{C, T}; rtol::Real=Base.rtoldefault(S, T), atol::Real=0, ztol::Real=1e-6)
    m = length(p.ps)
    if length(q.ps) != m
        false
    else
        permcomp((i, j) -> isapprox(p.ps[i], q.ps[j], rtol=rtol, atol=atol, ztol=ztol), m)
    end
end

isapprox{C, S, T, U, V}(p::RationalPoly{C, S, T}, q::RationalPoly{C, U, V}; rtol::Real=Base.rtoldefault(promote_op(*, U, T), promote_op(*, S, V)), atol::Real=0, ztol::Real=1e-6) = isapprox(p.num*q.den, q.num*p.den, rtol=rtol, atol=atol, ztol=ztol)
isapprox{C, S, T, U}(p::RationalPoly{C, S, T}, q::TermContainer{C, U}; rtol::Real=Base.rtoldefault(promote_op(*, U, T), S), atol::Real=0, ztol::Real=1e-6) = isapprox(p.num, q*p.den, rtol=rtol, atol=atol, ztol=ztol)
isapprox{C, S, T, U}(p::TermContainer{C, U}, q::RationalPoly{C, S, T}; rtol::Real=Base.rtoldefault(promote_op(*, U, T), S), atol::Real=0, ztol::Real=1e-6) = isapprox(p*q.den, q.num, rtol=rtol, atol=atol, ztol=ztol)
isapprox{C}(p::RationalPoly{C}, q; rtol::Real=Base.rtoldefault(promote_op(*, U, T), S), atol::Real=0, ztol::Real=1e-6) = isapprox(p, TermContainer{C}(q), rtol=rtol, atol=atol, ztol=ztol)
isapprox{C}(p, q::RationalPoly{C}; rtol::Real=Base.rtoldefault(promote_op(*, U, T), S), atol::Real=0, ztol::Real=1e-6) = isapprox(TermContainer{C}(p), q, rtol=rtol, atol=atol, ztol=ztol)

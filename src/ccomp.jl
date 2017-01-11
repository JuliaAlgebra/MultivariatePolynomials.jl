using Combinatorics
export isapproxzero

(==)(y, p::TermContainer) = TermContainer(y) == p
(==)(y::PolyType, p::TermContainer) = TermContainer(y) == p

function (==)(s::Term, t::Term)
    (s.α == t.α) && (iszero(s.α) || s.x == t.x)
end
function (==)(t::Term, p::VecPolynomial)
    if iszero(t.α)
        isempty(p.a)
    else
        length(p.a) == 1 && p.a[1] == t.α && p.x[1] == t.x
    end
end
(==)(p::VecPolynomial, t::Term) = t == p
function (==)(p::VecPolynomial, q::VecPolynomial)
    # terms should be sorted and without zeros
    for (tp,tq) in zip(p,q)
        if tp.x != tq.x
            # There should not be zero terms
            # FIXME if p is Term, it could be zero :/
            @assert tp.α != 0
            @assert tq.α != 0
            return false
        end
        if tp.α != tq.α
            return false
        end
    end
    true
end

(==)(p::RationalPoly, q::RationalPoly) = p.num*q.den == q.num*p.den
(==)(p, q::RationalPoly) = p*q.den == q.num

function isless(x::Vector, y::Vector)
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

# graded lex ordering
function isless(x::Monomial, y::Monomial)
    mycomp(x, y) < 0
end
isless(x::Monomial, y::PolyVar) = isless(x, Monomial(y))
isless(x::PolyVar, y::Monomial) = isless(Monomial(x), y)

function isapproxzero(x; ztol::Real=1e-6)
    -ztol < x < ztol
end

function isapproxzero(p::VecPolynomial; ztol::Real=1e-6)
    isapprox(p, zero(p), ztol=ztol)
end

function isapproxzero(p::RationalPoly; ztol::Real=1e-6)
    isapproxzero(p.num, ztol=ztol)
end

function isapprox{S,T}(p::VecPolynomial{S}, q::VecPolynomial{T}; rtol::Real=Base.rtoldefault(S, T), atol::Real=0, ztol::Real=1e-6)
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

function isapprox{S, T}(s::Term{S}, t::Term{T}; rtol::Real=Base.rtoldefault(S, T), atol::Real=0, ztol::Real=1e-6)
    s.x == t.x && isapprox(s.α, t.α, rtol=rtol, atol=atol)
end

function isapprox{S,T}(p::SOSDecomposition{S}, q::SOSDecomposition{T}; rtol::Real=Base.rtoldefault(S, T), atol::Real=0, ztol::Real=1e-6)
    m = length(p.ps)
    if length(q.ps) != m
        false
    else
        for σ in permutations(1:m)
            ok = true
            for i in 1:m
                if !isapprox(p.ps[i], q.ps[σ[i]], rtol=rtol, atol=atol, ztol=ztol)
                    ok = false
                    break
                end
            end
            if ok
                return true
            end
        end
        false
    end
end

isapprox{S,T,U,V}(p::RationalPoly{S,T}, q::RationalPoly{U,V}; rtol::Real=Base.rtoldefault(promote_op(*, U, T), promote_op(*, S, V)), atol::Real=0, ztol::Real=1e-6) = isapprox(p.num*q.den, q.num*p.den, rtol=rtol, atol=atol, ztol=ztol)
isapprox{S,T,U}(p::RationalPoly{S,T}, q::TermContainer{U}; rtol::Real=Base.rtoldefault(promote_op(*, U, T), S), atol::Real=0, ztol::Real=1e-6) = isapprox(p.num, q*p.den, rtol=rtol, atol=atol, ztol=ztol)
isapprox{S,T,U}(p::TermContainer{U}, q::RationalPoly{S,T}; rtol::Real=Base.rtoldefault(promote_op(*, U, T), S), atol::Real=0, ztol::Real=1e-6) = isapprox(p*q.den, q.num, rtol=rtol, atol=atol, ztol=ztol)
isapprox(p::RationalPoly, q; rtol::Real=Base.rtoldefault(promote_op(*, U, T), S), atol::Real=0, ztol::Real=1e-6) = isapprox(p, TermContainer(q), rtol=rtol, atol=atol, ztol=ztol)
isapprox(p, q::RationalPoly; rtol::Real=Base.rtoldefault(promote_op(*, U, T), S), atol::Real=0, ztol::Real=1e-6) = isapprox(TermContainer(p), q, rtol=rtol, atol=atol, ztol=ztol)

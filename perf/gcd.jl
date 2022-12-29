using BenchmarkTools

function bench0(x, T=Rational{BigInt})
    o = one(T)
    p = o * x^2 + 2o * x + o
    q = o * x + o
    @benchmark gcd($p, $q)
end

function bench1(x, y, z, T=Rational{BigInt})
    o = one(T)
    p = (o * x + o * y^2) * (o * z^3 + o * y^2 + o * x)
    q = (o * x + o * y + o * z) * (o * x^2 + o * y)
    @benchmark gcd($p, $q)
end

function bench2(x, y, z, t)
    p = x * y + 3 * (z * t)
    q = (p + 1) * p
    @benchmark gcd($p, $q)
end

function bench3(x, y, z, t)
    c1 = 10*(x * z + x)
    c2 = 2*(x^2 + z)
    c3 = 2*(2 - z  )
    c4 = 20*(x * z^2)
    e1 = 0
    e2 = 5
    e3 = 7
    e4 = 10
    p = c1 * y^e1 + c2 * y^e2 + c3 * y^e3 + c4 * y^e4
    q = prod(i->p + i, 0:3)
    @benchmark for i in 0:3
        gcd($p + i, $q)
    end
end

import SIMDPolynomials
function bench_SP()
    x, y, z, t = [SIMDPolynomials.PackedMonomial{4,7}(i) for i in 0:3]
    b0 = bench0(x)
    b1 = bench1(x, y, z)
    b2 = bench2(x, y, z, t)
    b3 = bench3(x, y, z, t)
    return b0, b1, b2, b3
end

import DynamicPolynomials
function bench_DP()
    DynamicPolynomials.@polyvar x y z t
    b0 = bench0(x)
    b1 = bench1(x, y, z)
    b2 = bench2(x, y, z, t)
    b3 = bench3(x, y, z, t)
    return b0, b1, b2, b3
end

import TypedPolynomials
function bench_TP()
    TypedPolynomials.@polyvar x y z t
    b0 = bench0(x)
    b1 = bench1(x, y, z)
    b2 = bench2(x, y, z, t)
    b3 = bench3(x, y, z, t)
    return b0, b1, b2, b3
end

include("table.jl")

function bench()
    bs0, bs1, bs2, bs3 = bench_SP()
    bd0, bd1, bd2, bd3 = bench_DP()
    bt0, bt1, bt2, bt3 = bench_TP()
    prettyprint(bs0, bd0, bt0)
    prettyprint(bs1, bd1, bt1)
    prettyprint(bs2, bd2, bt2)
    prettyprint(bs3, bd3, bt3)
end

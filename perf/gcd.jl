using BenchmarkTools

function bench0(x, T)
    o = one(T)
    p = o * x^2 + 2o * x + o
    q = o * x + o
    @benchmark gcd($p, $q)
end

function bench1(x, y, z, T)
    o = one(T)
    p = (o * x + o * y^2) * (o * z^3 + o * y^2 + o * x)
    q = (o * x + o * y + o * z) * (o * x^2 + o * y)
    @benchmark gcd($p, $q)
end

function bench2(x, y, z, t, T)
    p = x * y + T(3) * (z * t)
    q = (p + T(1)) * p
    @benchmark gcd($p, $q)
end

function bench3(x, y, z, t, T)
    c1 = T(10)*(x * z + x)
    c2 = T(2)*(x^2 + z)
    c3 = T(2)*(2 - z)
    c4 = T(20)*(x * z^2)
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
function bench_SP(T)
    x, y, z, t = [SIMDPolynomials.PackedMonomial{4,7}(i) for i in 0:3]
    b0 = bench0(x, T)
    b1 = bench1(x, y, z, T)
    b2 = bench2(x, y, z, t, T)
    b3 = bench3(x, y, z, t, T)
    return b0, b1, b2, b3
end

import DynamicPolynomials
function bench_DP(T)
    DynamicPolynomials.@polyvar x y z t
    b0 = bench0(x, T)
    b1 = bench1(x, y, z, T)
    b2 = bench2(x, y, z, t, T)
    b3 = bench3(x, y, z, t, T)
    return b0, b1, b2, b3
end

import TypedPolynomials
function bench_TP(T)
    TypedPolynomials.@polyvar x y z t
    b0 = bench0(x, T)
    b1 = bench1(x, y, z, T)
    b2 = bench2(x, y, z, t, T)
    b3 = bench3(x, y, z, t, T)
    return b0, b1, b2, b3
end

include("table.jl")

function bench(T)
    bs0, bs1, bs2, bs3 = bench_SP(T)
    bd0, bd1, bd2, bd3 = bench_DP(T)
    bt0, bt1, bt2, bt3 = bench_TP(T)
    println("### Benchmark 0")
    println()
    prettyprint(bs0, bd0, bt0)
    println()
    println("### Benchmark 1")
    println()
    prettyprint(bs1, bd1, bt1)
    println()
    println("### Benchmark 2")
    println()
    prettyprint(bs2, bd2, bt2)
    println()
    println("### Benchmark 3")
    println()
    prettyprint(bs3, bd3, bt3)
    println()
end

function bench()
    for T in [Int, BigInt, Rational{BigInt}]
        println("## $T")
        println()
        bench(T)
    end
end

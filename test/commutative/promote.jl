using Test
import MutableArithmetics as MA
using MultivariatePolynomials
import MultivariatePolynomials as MP
using LinearAlgebra

@testset "Promotion" begin
    Mod.@polyvar x y
    @inferred x * y + x
    @test [x, x * y + x, x] isa Vector{<:AbstractPolynomial{Int}}
    @test eltype([1, x / y, x]) <:
          RationalPoly{<:AbstractTerm{Int},<:AbstractTerm{Int}}
    @test [(x^2-2x+2) x; x x^2] isa Matrix{<:AbstractPolynomial{Int}}
    @test [2.0x, 3x] isa Vector{<:AbstractTerm{Float64}}
    @inferred Any[x*y, x+y]
    @test Any[x*y, x+y] isa Vector{Any}
    @test [x * y, x + y] isa Vector{<:AbstractPolynomial{Int}}
    @test [2x * y, x + y] isa Vector{<:AbstractPolynomial{Int}}
    @test eltype([2.0x, x / y]) <:
          RationalPoly{<:AbstractTerm{Float64},<:AbstractTerm{Int}}
    @test eltype([2.0x, x / y, 1y]) <:
          RationalPoly{<:AbstractTerm{Float64},<:AbstractTerm{Int}}
    @test eltype([2x + y, x / 2.0y, x + 1y]) <:
          RationalPoly{<:AbstractPolynomial{Int},<:AbstractTerm{Float64}}
    @test eltype([-1 / x^2, 1]) <:
          RationalPoly{<:AbstractTerm{Int},<:AbstractTerm{Int}}

    X = [x, y]
    Y = [1 2; 3 4] * X
    @test Y[1] == x + 2y
    @test Y[2] == 3x + 4y
    Y = X' * [1 2; 3 4]
    @test Y[1] == x + 3y
    @test Y[2] == 2x + 4y
    @test dot(X, [1 2; 3 4] * X) == x^2 + 5x * y + 4y^2

    function _t(a, b, T)
        if VERSION < v"1.6"
            @test typeof(@inferred vcat(a, b)) in [Vector{T}, Vector{Any}]
            @test typeof(@inferred vcat(b, a)) in [Vector{T}, Vector{Any}]
        else
            @test typeof(@inferred vcat(a, b)) == Vector{T}
            @test typeof(@inferred vcat(b, a)) == Vector{T}
        end
    end

    function __test(a, PT, TT, MT)
        @test typeof(a) == Vector{MT(Int)}
        _t(a, 1, Any)
        _t(a, x, MT(Int))
        _t(a, x^2, MT(Int))
        _t(a, 1x, TT(Int))
        _t(a, 1.0x, TT())
        _t(a, x + y, PT(Int))
        return _t(a, 1.0x + y, PT())
    end

    function _test(a, af, TT)
        __test(a, apl, TT, TT)
        # `x isa _APL{Int}` so here we don't have `Float64`:
        pt() = apl()
        pt(T) = apl()
        tt() = TT()
        tt(T) = TT()
        __test(af, pt, tt, tt)
        @test typeof(af) == Vector{TT()}
        return _t(af, a, TT())
    end

    apl() = MP._APL
    apl(T::Type) = MP._APL{T}
    p = [i == 1 ? x + y : x for i in 1:2]
    pf = [i == 1 ? 1.0x + y : x for i in 1:2]
    _test(p, pf, apl)

    atl() = MP.AbstractTermLike
    atl(T::Type) = MP.AbstractTermLike{T}
    t = [i == 1 ? 2x : x for i in 1:2]
    tf = [i == 1 ? 2.0x : x for i in 1:2]
    _test(t, tf, atl)
    _t(p, t, apl(Int))
    _t(pf, t, apl())
    _t(p, tf, apl())
    _t(pf, tf, apl())

    __pt = Base.promote_typejoin(typeof(x + 2), typeof(x + 2.0))
    _pt() = __pt
    _pt(::Type) = __pt
    __test([i == 1 ? x + 2 : x + 2.0 for i in 1:2], _pt, _pt, _pt)

    __tt = Base.promote_typejoin(typeof(2x), typeof(2.0x))
    _tt() = __tt
    _tt(::Type) = __tt
    __test([i == 1 ? 2x : 2.0x for i in 1:2], _pt, _tt, _tt)

    aml(args...) = MP.AbstractMonomialLike
    a = [i == 1 ? x^2 : x for i in 1:2]
    __test(a, apl, atl, aml)
    _t(a, t, atl(Int))
    _t(a, tf, atl())
    _t(a, p, apl(Int))
    _t(a, pf, apl())
end

struct A end
Base.zero(::Type{A}) = A()
Base.:*(::A, ::A) = B()
struct B end
Base.zero(::Type{B}) = B()
Base.:+(::B, ::B) = C()
struct C end
Base.:+(::B, ::C) = C()
Base.:+(::C, ::B) = C()
Base.:+(::C, ::C) = C()
@testset "promote_operation with tricky types" begin
    Mod.@polyvar x
    @test MA.promote_operation(
        *,
        polynomial_type(x, A),
        polynomial_type(x, A),
    ) == polynomial_type(x, C)
    @test MA.promote_operation(*, polynomial_type(x, A), A) ==
          polynomial_type(x, B)
    @test MA.promote_operation(*, A, polynomial_type(x, A)) ==
          polynomial_type(x, B)
end

function __promote_prod(::Type{A}, ::Type{B}, ::Type{C}) where {A,B,C}
    @test MA.promote_operation(*, A, B) == C
    @test MA.promote_operation(*, B, A) == C
end

@testset "promote_operation with Rational" begin
    Mod.@polyvar x
    V = typeof(x)
    M = monomial_type(V)
    T = term_type(V, Int)
    P = polynomial_type(V, Float64)
    function _promote_prod(::Type{A}, ::Type{B}, ::Type{C}) where {A,B,C}
        __promote_prod(A, B, C)
        __promote_prod(RationalPoly{A,B}, RationalPoly{B,A}, RationalPoly{C,C})
        __promote_prod(RationalPoly{A,A}, RationalPoly{B,B}, RationalPoly{C,C})
        for U in [V, M, T, P]
            __promote_prod(A, RationalPoly{B,U}, RationalPoly{C,U})
        end
    end
    _promote_prod(V, V, M)
    for U in [V, M]
        _promote_prod(U, M, M)
    end
    for U in [V, M, T]
        _promote_prod(U, T, T)
    end
    for U in [V, M, T, P]
        _promote_prod(U, P, P)
    end
end

@testset "promote_operation with polynomial coefficient" begin
    Mod.@polyvar x
    Mod.@polyvar y
    X = typeof(x)
    MX = monomial_type(X)
    TX = term_type(X, Int)
    PX = polynomial_type(X, Int)
    Y = typeof(y)
    TXY = term_type(Y, PX)
    PXY = polynomial_type(Y, PX)
    TY = term_type(Y, Int)
    PY = polynomial_type(Y, Int)
    for T in [X, MX, TX, PX]
        __promote_prod(TY, TY, TY)
        __promote_prod(TXY, TY, TXY)
        __promote_prod(TXY, TXY, TXY)
        __promote_prod(PY, PY, PY)
        __promote_prod(PXY, PY, PXY)
        __promote_prod(PXY, PXY, PXY)
    end
end

import StarAlgebras as SA

@testset "ExponentMap" begin
    # Vector exponent map
    map_v = MP.ExponentMap([2, 3], 4)
    @test map_v([5, 7]) == [0, 5, 7, 0]

    # Tuple exponent map
    map_t = MP.ExponentMap((2, 3), Val(4))
    @test map_t((5, 7)) == (0, 5, 7, 0)

    # Identity-like map (all slots filled)
    map_id = MP.ExponentMap([1, 2], 2)
    @test map_id([3, 4]) == [3, 4]
end

@testset "promote_variables_with_maps" begin
    # Same variables: no maps needed
    Mod.@polyvar a b
    vars_ab = MP.variables(a + b)
    (va, ma), (vb, mb) = MP.promote_variables_with_maps(vars_ab, vars_ab)
    @test va == vars_ab
    @test vb == vars_ab
    @test ma === nothing
    @test mb === nothing

    # Different variables from separate @polyvar calls: maps should be created
    Mod.@polyvar x
    Mod.@polyvar y
    vars_x = MP.variables(x^2)
    vars_y = MP.variables(y^2)
    @test vars_x != vars_y
    (vx, mx), (vy, my) = MP.promote_variables_with_maps(vars_x, vars_y)
    @test length(vx) == 2
    @test vx == vy  # both should be the same merged variable set
    @test mx isa MP.ExponentMap
    @test my isa MP.ExponentMap

    exp_x = MP.exponents(x^3)
    new_exp_x = mx(exp_x)
    @test sum(new_exp_x) == sum(exp_x)
    @test length(new_exp_x) == length(vx)

    exp_y = MP.exponents(y^2)
    new_exp_y = my(exp_y)
    @test sum(new_exp_y) == sum(exp_y)
    @test length(new_exp_y) == length(vy)

    # The expanded exponents should not overlap (different slots)
    @test all(new_exp_x .* new_exp_y .== 0)
end

@testset "promote_bases_with_maps" begin
    # Same variables: no promotion needed (maybe_promote is a no-op)
    Mod.@polyvar a b
    p = a^2 + a
    q = a^2 + 2a
    (pp, mp), (qq, mq) = SA.promote_bases_with_maps(p, q)
    @test pp === p
    @test qq === q
    @test mp === nothing
    @test mq === nothing

    # Same variables via promote_bases
    pp2, qq2 = SA.promote_bases(p, q)
    @test pp2 === p
    @test qq2 === q
end

@testset "promote_bases end-to-end" begin
    # Monomials with different variables
    Mod.@polyvar x
    Mod.@polyvar y
    px, py = SA.promote_bases(x^2, y^2)
    @test MP.variables(px) == MP.variables(py)
    @test length(MP.variables(px)) == 2
    @test MP.degree(px) == 2
    @test MP.degree(py) == 2
    # The original exponent structure is preserved
    @test MP.degree(px, x) == 2
    @test MP.degree(px, y) == 0
    @test MP.degree(py, x) == 0
    @test MP.degree(py, y) == 2

    # Variables with different variables
    px2, py2 = SA.promote_bases(x, y^3)
    @test MP.variables(px2) == MP.variables(py2)
    @test length(MP.variables(px2)) == 2
    @test MP.degree(px2) == 1
    @test MP.degree(py2) == 3

    # promote_bases_with_maps returns maps for different variables
    (rx, mx), (ry, my) = SA.promote_bases_with_maps(x^2, y^2)
    @test MP.variables(rx) == MP.variables(ry)
    @test mx isa MP.ExponentMap
    @test my isa MP.ExponentMap

    # Term with Monomial
    t = 3x^2
    m = y^2
    pt, pm = SA.promote_bases(t, m)
    @test MP.variables(pt) == MP.variables(pm)
    @test length(MP.variables(pt)) == 2
    @test MP.coefficient(pt) == 3
    @test MP.degree(pt) == 2
    @test MP.degree(pt, x) == 2
    @test MP.degree(pt, y) == 0
    @test MP.degree(pm, x) == 0
    @test MP.degree(pm, y) == 2

    # Term with Term
    t1 = 2x^2
    t2 = 5y^3
    pt1, pt2 = SA.promote_bases(t1, t2)
    @test MP.variables(pt1) == MP.variables(pt2)
    @test MP.coefficient(pt1) == 2
    @test MP.coefficient(pt2) == 5
    @test MP.degree(pt1, x) == 2
    @test MP.degree(pt1, y) == 0
    @test MP.degree(pt2, x) == 0
    @test MP.degree(pt2, y) == 3

    # Term with same variables (no promotion needed)
    t3 = 4x^2
    t4 = 7x^3
    pt3, pt4 = SA.promote_bases(t3, t4)
    @test MP.coefficient(pt3) == 4
    @test MP.coefficient(pt4) == 7

    # Polynomial with Monomial
    p = 2x^2 + 3x
    m = y^2
    pp, pm = SA.promote_bases(p, m)
    @test MP.variables(pp) == MP.variables(pm)
    @test length(MP.variables(pp)) == 2
    @test MP.degree(pm, y) == 2
    @test MP.degree(pm, x) == 0

    # Polynomial with Polynomial (different variables)
    p1 = x^2 + 2x
    p2 = y^3 + 3y
    pp1, pp2 = SA.promote_bases(p1, p2)
    @test MP.variables(pp1) == MP.variables(pp2)
    @test length(MP.variables(pp1)) == 2
end

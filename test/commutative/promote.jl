if !@isdefined(Mod)
    include("../SimplePolynomials.jl")
    import .SimplePolynomials as Mod
    include("../utils.jl")
end

using Test
import MutableArithmetics as MA
using MultivariatePolynomials
import MultivariatePolynomials as MP
import StarAlgebras as SA
using LinearAlgebra

@testset "Fused promotion preserves operand bases" begin
    Mod.@polyvar x y z
    f = x + 1
    for op in (MA.add_mul, MA.sub_mul), a in (2, 2y, y + 1), b in (3z, z + 1)
        expected = op(f, a, b)
        @test (@inferred MA.promote_operation(
            op,
            typeof(f),
            typeof(a),
            typeof(b),
        )) === typeof(expected)
        @test (@inferred MA.promote_operation(
            op,
            typeof(f),
            typeof(b),
            typeof(a),
        )) === typeof(op(f, b, a))
    end
end

@testset "Promotion" begin
    Mod.@polyvar x y
    @testset "Promotion with unrelated types" begin
        p = x + y
        q = x + 1.0y
        @test (@inferred promote_type(typeof(p), typeof(q))) == typeof(q)
        for value in (:(ztol = 1e-7), :symbol, "string")
            for poly in (x, x^2, 2x, p, x / y, p / y)
                @test promote_type(typeof(value), typeof(poly)) == Any
                @test promote_type(typeof(poly), typeof(value)) == Any
                @test [value, poly] isa Vector{Any}
                @test [poly, value] isa Vector{Any}
            end
            # Julia 1.13's `@test isapprox(p, q, ztol = ...)` builds
            # `[kwargs_expr, p, q]` for diagnostics, even when the test passes.
            # This promotes the keyword Expr with the polynomial-like arguments
            # and must not turn unrelated types into polynomial coefficients.
            # Julia 1.10 promotes p and q first; newer Julia versions start
            # with value and p. Both orders must give Any.
            @test Base.promote_typeof(value, p, q) == Any
            @test [value, p, q] isa Vector{Any}
        end
    end
    @inferred x * y + x
    @test [x, x * y + x, x] isa Vector{<:AbstractPolynomial{Int}}
    @test eltype([1, x / y, x]) <:
          RationalPoly{<:AbstractTerm{Int},<:AbstractTerm{Bool}}
    @test [(x^2-2x+2) x; x x^2] isa Matrix{<:AbstractPolynomial{Int}}
    @test [2.0x, 3x] isa Vector{<:AbstractTerm{Float64}}
    @inferred Any[x*y, x+y]
    @test Any[x*y, x+y] isa Vector{Any}
    @test [x * y, x + y] isa Vector{<:AbstractPolynomial{Int}}
    @test [2x * y, x + y] isa Vector{<:AbstractPolynomial{Int}}
    @test eltype([2.0x, x / y]) <:
          RationalPoly{<:AbstractTerm{Float64},<:AbstractTerm{Bool}}
    @test eltype([2.0x, x / y, 1y]) <:
          RationalPoly{<:AbstractTerm{Float64},<:AbstractTerm{Bool}}
    @test eltype([2x + y, x / 2.0y, x + 1y]) <:
          RationalPoly{<:AbstractPolynomial{Int},<:AbstractTerm{Float64}}
    @test eltype([-1 / x^2, 1]) <:
          RationalPoly{<:AbstractTerm{Int},<:AbstractTerm{Bool}}

    X = [x, y]
    Y = [1 2; 3 4] * X
    @test Y[1] == x + 2y
    @test Y[2] == 3x + 4y
    Y = X' * [1 2; 3 4]
    @test Y[1] == x + 3y
    @test Y[2] == 2x + 4y
    @test dot(X, [1 2; 3 4] * X) == x^2 + 5x * y + 4y^2

    # Concatenation promotes concrete terms and polynomial coefficients.
    for (a, b, T) in (
        ([x^2], x, monomial_type(x)),
        ([x^2], 2x, term_type(x, Int)),
        ([2x], 2.0x, term_type(x, Float64)),
        ([2x, 2.0x], 3x, term_type(x, Float64)),
        ([x + y], 2.0x, polynomial_type(x * y, Float64)),
        ([x + 2, x + 2.0], y, polynomial_type(x * y, Float64)),
    )
        ab = @inferred vcat(a, b)
        ba = @inferred vcat(b, a)
        @test ab isa Vector{T}
        @test ba isa Vector{T}
        @test ab[1:(end-1)] == ba[2:end] == a
        @test ab[end] == ba[1] == b
    end

    # Concrete type promote_type tests
    @test promote_type(typeof(x), typeof(y)) <: AbstractMonomialLike
    @test promote_type(typeof(2x), typeof(3y)) <: AbstractTerm{Int}
    @test promote_type(typeof(2x), typeof(3.0y)) <: AbstractTerm{Float64}
    @test promote_type(typeof(x + y), typeof(x)) <: AbstractPolynomialLike{Int}
    @test promote_type(typeof(2x + y), typeof(3.0y)) <:
          AbstractPolynomialLike{Float64}

    for p in (x, x^2, 2x, x + y)
        @test promote_type(typeof(p), String) === Any
    end
end

struct A end
Base.zero(::Type{A}) = A()
Base.:*(::A, ::A) = B()
struct B end
Base.zero(::Type{B}) = B()
# Algebra products multiply coefficient products by integer basis weights.
Base.:*(::Int, ::B) = B()
Base.:+(::B, ::B) = C()
struct C end
Base.:+(::B, ::C) = C()
Base.:+(::C, ::B) = C()
Base.:+(::C, ::C) = C()
# `A` is not a scalar by default, register it for polynomial coefficients
SA.@scalar_arithmetic A
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
    a = x + 1
    p = polynomial([a, x], [y, one(y)])
    expected = polynomial([a * a, a * x], [y, one(y)])
    @test (@inferred a * p) == expected
    @test (@inferred p * a) == expected
    @test coefficient_type(a * p) === typeof(a)
    for (left, right) in ((a, p), (p, a))
        output = zero(p)
        @test (@inferred MA.operate_to!!(output, *, left, right)) === output
        @test output == expected
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

    for (a, b) in ((x, y), (x, y^2), (2x, y^2), (x + 1, y^2))
        aa, bb = @inferred SA.promote_bases(a, b)
        rb, ra = @inferred SA.promote_bases(b, a)
        @test variables(aa) == variables(bb) == variables(ra) == variables(rb)
        @test aa == ra == a
        @test bb == rb == b
    end
end

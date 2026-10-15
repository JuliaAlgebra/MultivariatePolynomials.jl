if !@isdefined(Mod)
    import DynamicPolynomials as Mod
    include("../utils.jl")
end

import StarAlgebras as SA

@testset "Aligned exponents" begin
    Mod.@polyvar x y
    p = x^2 * y
    e = @inferred exponents(p, variables(p))
    @test e == exponents(p)
    @test typeof(e) == typeof(exponents(p))
    if e isa AbstractVector
        @test e !== exponents(p)
    end
    @test collect(@inferred exponents(x^2, variables(p))) == [2, 0]
    @test_throws ArgumentError exponents(p, variables(x))
end

@testset "Numeric monomial conversion" begin
    Mod.@polyvar x
    M = monomial_type(x)
    for value in (1, 1.0, 1 // 1)
        @test (@inferred convert(M, value)) isa M
        @test convert(M, value) == one(x)
    end
    for value in (0, 2, 0.5)
        @test_throws InexactError convert(M, value)
    end
end

function effective_variables_test()
    Mod.@polyvar x
    Mod.@polyvar y[1:7]
    T = variable_union_type(x)
    @test x isa T
    @test y[2] isa T
    @test T[x, y[2]] == @inferred effective_variables(x * y[2])
    @test T[x, y[2]] == @inferred effective_variables(y[2] * x)
    @test T[x] == @inferred effective_variables(x * y[2]^0)
    @test T[x] == @inferred effective_variables(y[2]^0 * x)
    @test T[y[2]] == @inferred effective_variables(x^0 * y[2])
    @test T[y[2]] == @inferred effective_variables(y[2] * x^0)
    @test T[x, y[2]] == @inferred effective_variables(y[3]^0 * x * y[2])
    @test T[x, y[2]] == @inferred effective_variables(y[2] * y[3]^0 * x)
    @test T[x, y[2]] == @inferred effective_variables(y[2] * x * y[3]^0)
end

function map_exponents_test()
    Mod.@polyvar x y
    a = x^2
    b = x * y
    shared = one(a)
    @test_throws ArgumentError MP.map_exponents!(+, a, b)
    @test a == x^2
    @test collect(variables(shared)) == [x]

    a, b = SA.promote_bases(a, b)
    c = @inferred MP.map_exponents!(+, a, b)
    @test (c === a) == (exponents(a) isa AbstractVector)
    @test c == x^3 * y
    @test variables(c) == variables(b)
    @test b == x * y

    a = x^3
    @test_throws ArgumentError MP.map_exponents_to!(a, -, b, b)
    @test a == x^3
    a, b = SA.promote_bases(a, b)
    d = @inferred MP.map_exponents_to!(a, -, b, b)
    @test (d === a) == (exponents(a) isa AbstractVector)
    @test isone(d)
    @test variables(d) == variables(b)

    # Smaller operand bases are aligned without changing the destination.
    output = one(x * y)
    vars = variables(output)
    output = @inferred MP.map_exponents_to!(output, +, x^2, y)
    @test output == x^2 * y
    @test variables(output) === vars
    output = @inferred MP.map_exponents!(+, output, y)
    @test output == x^2 * y^2

    # Either input may alias the destination.
    a, b = x^2 * y, x * y^3
    b = @inferred MP.map_exponents_to!(b, +, a, b)
    @test b == x^3 * y^4
    @test a == x^2 * y
    a = @inferred MP.map_exponents_to!(a, -, a, a)
    @test isone(a)
    @test (@inferred MP.map_exponents!(div, x^4 * y^6, x^2 * y^3)) == x^2 * y^2

    let output = output, a = x^2 * y, b = x * y^3
        alloc_test(() -> MP.map_exponents_to!(output, +, a, b), 0)
        alloc_test(() -> MP.map_exponents!(+, output, a), 0)
    end
    return
end

@testset "Monomial" begin
    Mod.@polyvar x

    @test zero(x^2) == 0
    @test zero(x^2) isa AbstractPolynomial{Int}
    @inferred zero(x^2)
    @test (@inferred one(x^2)) == 1
    @test one(x^2) isa AbstractMonomial
    @test (@inferred one(typeof(x^2))) == 1
    @test one(typeof(x^2)) isa AbstractMonomial
    # @test ordering(x) === GradedLex #FIXME add it back once it's implemented

    Mod.@polyvar y[1:7]
    mono = y[1] * y[3] * y[5] * y[7]
    alloc_test(() -> convert(typeof(mono), mono), 0)
    @test issorted(variables(y[2] * mono), rev = true)
    @test issorted(variables(mono * y[4]), rev = true)
    @test issorted(variables(y[6] * mono), rev = true)

    @test nvariables(monomial_vector([x^2, prod(y[2:4])])) == 4

    @test nterms(x^2) == 1
    @test @inferred(terms(x^2)) == [x^2]

    @test degree(x * y[2]^2, x) == 1
    @test degree(x * y[2]^2, y[1]) == 0
    @test degree(x * y[2]^2, y[2]) == 2

    @test_throws InexactError variable(x^2)
    @test_throws InexactError variable(x * y[1])
    @test_throws InexactError variable(constant_monomial(typeof(x)))

    @test x != constant_monomial(typeof(x))
    @test constant_monomial(typeof(x)) != x

    m = x^2
    @test x != m
    @test m != x
    typetests(m)
    typetests([x^2, x^3])
    @test (@inferred polynomial(m)) isa AbstractPolynomial{Int}
    @test (@inferred polynomial(m, Float64)) isa AbstractPolynomial{Float64}

    @test variable(x^1) == x
    @test variable(x^1) isa AbstractVariable
    @test variable(x^2 + x - x^2) == x
    @test variable(x^2 + x - x^2) isa AbstractVariable
    @test variable(1.0x) == x
    @test variable(1.0x) isa AbstractVariable
    @test_throws InexactError variable(x + 2x) == x

    @test monic(x^2) == x^2

    @test MP.div_multiple(2x^2 * y[1]^3, x * y[1]^2) == 2x * y[1]

    @test transpose(x) == x
    @test adjoint(x) == x
    @test transpose(x^2) == x^2
    @test adjoint(x^2) == x^2

    @testset "Effective variables" begin
        effective_variables_test()
    end
    @testset "map_exponents" begin
        map_exponents_test()
    end
end

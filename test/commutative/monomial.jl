if !@isdefined(Mod)
    include("../SimplePolynomials.jl")
    import .SimplePolynomials as Mod
    include("../utils.jl")
end

import StarAlgebras as SA
import MutableArithmetics as MA

@testset "Monomial bases" begin
    Mod.@polyvar x y
    full = @inferred MP.FullBasis{MP.Monomial}(x * y)
    basis = MP.SubBasis{MP.Monomial}([1, x, y])
    @test basis == full[monomial_vector([y, 1, x])]
    @test collect(basis) == monomial_vector([1, y, x])
    @test MP.keys_as_monomials(basis) == monomial_vector([1, y, x])
    @test MP.explicit_basis_covering(full, basis) == basis
    @test MP.explicit_basis_type(typeof(full)) === typeof(basis)
    @test isempty(MP.empty_basis(typeof(basis)))
    @test MP.maxdegree_basis(full, 2) ==
          MP.SubBasis{MP.Monomial}(monomials(variables(x * y), 0:2))

    vars = SA.object(MP.algebra(full))
    other_vars = MP.Variables{MP.Monomial}(MA.copy_if_mutable(variables(vars)))
    @test vars == other_vars
    @test hash(vars) == hash(other_vars)
    @test one(vars) == one(x * y)
    @test MP.BasisElement{MP.Monomial}(x * y) == x * y
    @test SA.star(full, full[x*y]) == full[x*y]

    for b in (full, basis)
        @test variables(b) == variables(x * y)
        @test nvariables(b) == 2
        @test MP.variable_index(b, x) == 1
        @test MP.variable_index(b, y) == 2
        @test monomial_type(typeof(b)) === monomial_type(x * y)
        @test (@inferred MP.ordering(b)) === MP.Graded(MP.LexOrder())
        @test (@inferred MP.ordering(typeof(b))) === MP.ordering(x)
        @test polynomial_type(b, Int) === polynomial_type(x * y, Int)
        @test polynomial_type(typeof(b), Float64) ===
              polynomial_type(x * y, Float64)
        @test typeof(MP.algebra(b)) ===
              MA.promote_operation(MP.algebra, typeof(b))
        @test polynomial_type(typeof(MP.algebra(b)), Int) ===
              polynomial_type(x * y, Int)
        for c in (2, 0.5, Float64)
            a = @inferred MP.constant_algebra_element(b, c)
            value = c isa Type ? one(c) : c
            @test a == value
            @test coefficient_type(a) === typeof(value)
            @test typeof(a) ===
                  MP.constant_algebra_element_type(typeof(b), typeof(value))
        end
    end

    @testset "Finite and full coefficients" begin
        p = SA.algebra_element(i -> i^2, basis)
        @test SA.coeffs(p) == [1, 4, 9]
        @test MP.explicit_basis(p) === basis
        @test typeof(p) === MP.algebra_element_type(Vector{Int}, typeof(basis))
        q = @inferred MP.implicit(p)
        @test q == 1 + 4y + 9x
        @test SA.basis(q) === parent(basis)
        @test typeof(q) === MA.promote_operation(MP.implicit, typeof(p))
        @test MP.implicit(q) === q
        @test MP.explicit_basis(q) == basis
        @test polynomial(i -> i^2, basis) == q
        @test coefficients(q, basis) == [1, 4, 9]
        @test coefficients(q, full) === q
        @test MP.explicit_basis(2x) == MP.SubBasis{MP.Monomial}([x])
        @test MP.explicit_basis(x * y) == MP.SubBasis{MP.Monomial}([x * y])
        for a in (x, x * y, 2x)
            @test typeof(SA.algebra_element(a)) ===
                  MA.promote_operation(SA.algebra_element, typeof(a))
        end
        for a in (basis[3], 2basis[3], p, q)
            converted = MP.convert_basis(full, a)
            @test converted == a
            @test SA.basis(converted) === full
        end
        @test SA.coeffs(MP.convert_basis(basis, q)) == [1, 4, 9]

        linear = MP.SubBasis{MP.Monomial}([x, y])
        @test SA.coeffs([3, 7], linear, basis) == [0, 3, 7]
        @test SA.adjoint_coeffs([1, 4, 9], linear, basis) == [4, 9]
        @test SA.adjoint_coeffs!([99, 99], [1, 4, 9], linear, basis) == [4, 9]
        @test polynomial([1 2; 2 3], linear, Int) == y^2 + 4x * y + 3x^2
        empty = MP.empty_basis(typeof(linear))
        @test iszero(polynomial(zeros(Int, 0, 0), empty, Int))
    end

    @testset "Multiplicative structure" begin
        m = MP.algebra(full).mstructure
        xy = m(basis[3], basis[2], monomial_type(x * y))
        @test collect(SA.keys(xy)) == [x * y]
        @test collect(SA.values(xy)) == [1]

        p = x - x^2
        promoted, promoted_full = SA.promote_bases(p, full)
        @test promoted == p
        @test SA.basis(promoted) == full
        @test promoted_full === full
        @test variables(p) == variables(x)
        # Finite bases keep their keys while their variable sets are promoted.
        a = SA.algebra_element([2], MP.SubBasis{MP.Monomial}([x]))
        b = SA.algebra_element([3], MP.SubBasis{MP.Monomial}([y]))
        promoted_a, promoted_b = SA.promote_bases(a, b)
        @test promoted_a == a == 2x
        @test promoted_b == b == 3y
        @test variables(promoted_a) == variables(promoted_b) == variables(full)
        @test variables(a) == variables(x)
        @test variables(b) == variables(y)
    end

    @testset "Printing" begin
        for (b, expected) in (
            (full, "FullBasis{Monomial}([x, y])"),
            (basis, "SubBasis{Monomial}([1, y, x])"),
            (vars, "Monomial polynomials in the variables [x, y]"),
        )
            @test sprint(show, b) == expected
            @test sprint(print, b) == expected
            for mime in (MIME"text/plain"(), MIME"text/print"())
                @test sprint(show, mime, b) == expected
            end
        end
        quadratic = MP.SubBasis{MP.Monomial}([x^2, x * y, y^2])
        @test sprint(show, quadratic) == "SubBasis{Monomial}([y², xy, x²])"
        @test sprint(print, quadratic) == "SubBasis{Monomial}([y^2, x*y, x^2])"
    end
end

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
    # See https://github.com/JuliaAlgebra/DynamicPolynomials.jl/pull/105
    e = x * y^0
    f = x * y^2
    @test variables(e) == variables(f)
    g = @inferred MP.map_exponents!(div, e, f)
    if exponents(e) isa AbstractVector
        @test g === e
    end
    @test g == x
    @test variables(g) == variables(f)
    @test degree(g, x) == 1
    @test degree(g, y) == 0

    let output = output, a = x^2 * y, b = x * y^3
        alloc_test(() -> MP.map_exponents_to!(output, +, a, b), 0)
        alloc_test(() -> MP.map_exponents!(+, output, a), 0)
    end
    return
end

@testset "Monomial" begin
    Mod.@polyvar x

    @test zero(x^2) == 0
    @test zero(x^2) isa AbstractPolynomial{Bool}
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
    typetests(m, Bool)
    typetests([x^2, x^3], Bool)
    @test (@inferred polynomial(m)) isa AbstractPolynomial{Bool}
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

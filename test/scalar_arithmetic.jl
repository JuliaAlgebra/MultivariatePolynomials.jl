module TestScalarArithmetic

using Test
import DynamicPolynomials as DP
import TypedPolynomials as TP
import MultivariatePolynomials as MP
import StarAlgebras as SA
import MutableArithmetics as MA
import JuMP
import MathOptInterface as MOI
import LinearAlgebra

struct UnregisteredMutable <: MA.AbstractMutable end

function test_equal(p, q)
    @test MP.monomials(p) == MP.monomials(q)
    @test length(MP.coefficients(p)) == length(MP.coefficients(q))
    @test all(MA.isequal_canonical.(MP.coefficients(p), MP.coefficients(q)))
    return
end

function test_monomial(a, m)
    for product in ((@inferred a * m), (@inferred m * a))
        @test product isa MP.AbstractTerm{typeof(a)}
        @test MP.coefficient(product) === a
        @test MP.monomial(product) == m
    end
    @test MA.promote_operation(*, typeof(a), typeof(m)) === typeof(a * m)
    @test MA.promote_operation(*, typeof(m), typeof(a)) === typeof(m * a)
    for op in (+, -), left in (false, true)
        p = left ? (@inferred op(a, m)) : (@inferred op(m, a))
        A, B = left ? (typeof(a), typeof(m)) : (typeof(m), typeof(a))
        @test (@inferred MA.promote_operation(op, A, B)) === typeof(p)
        @test MA.isequal_canonical(
            MP.coefficient(p, m),
            convert(MP.coefficient_type(p), left ? op(1) : 1),
        )
        @test MA.isequal_canonical(
            MP.coefficient(p, one(m)),
            left ? op(a, 0) : op(0, a),
        )
    end
    for (op, add) in ((MA.add_mul, +), (MA.sub_mul, -)), left in (false, true)
        p = MP.polynomial(m + 2)
        expected = add(p, left ? a * (1m) : (1m) * a)
        output = convert(typeof(expected), MA.mutable_copy(p))
        result =
            left ? MA.operate!(op, output, a, m) : MA.operate!(op, output, m, a)
        @test result === output
        test_equal(result, expected)
        test_equal(p, m + 2)
    end
    test_equal(LinearAlgebra.dot(a, m), adjoint(a) * m)
    test_equal(LinearAlgebra.dot(m, a), adjoint(m) * a)
    return
end

function test_dots(a, monos)
    scalars = [a, a]
    for left in (false, true)
        x, y = left ? (scalars, monos) : (monos, scalars)
        expected = adjoint(x[1]) * y[1] + adjoint(x[2]) * y[2]
        test_equal(@inferred(LinearAlgebra.dot(x, y)), expected)
        test_equal(@inferred(MA.operate(LinearAlgebra.dot, x, y)), expected)
        @test_throws DimensionMismatch LinearAlgebra.dot(x, y[1:1])
    end
    return
end

scalar_array(p, a) = [p, a]

function test_promotion(a, p)
    P = @inferred promote_type(typeof(p), typeof(a))
    @test isconcretetype(P)
    @test (@inferred promote_type(typeof(a), typeof(p))) === P
    for entries in
        ((@inferred scalar_array(p, a)), reverse(@inferred scalar_array(a, p)))
        @test eltype(entries) === P
        original, constant = entries
        if p isa MP.RationalPoly
            C = MP.coefficient_type(numerator(original))
            test_equal(
                numerator(original),
                MP.map_coefficients(Base.Fix1(convert, C), numerator(p)),
            )
            test_equal(denominator(original), denominator(p))
            @test isone(denominator(constant))
            constant = numerator(constant)
        else
            C = MP.coefficient_type(original)
            test_equal(original, MP.map_coefficients(Base.Fix1(convert, C), p))
        end
        C = MP.coefficient_type(constant)
        test_equal(constant, MP.constant_term(convert(C, a), constant))
        @test MA.isequal_canonical(
            @inferred(convert(C, constant)),
            convert(C, a),
        )
    end
    return
end

function test_scalar(a, x)
    test_monomial(a, x)
    test_monomial(a, x^2)
    C = MA.promote_operation(+, Int, typeof(a))
    M = MP.monomial_type(x)
    @test (@inferred convert(M, one(C))) == one(x)
    @test_throws InexactError convert(M, zero(C))
    for p in (x, x^2, 2x, 2x + 3, x / (x + 1))
        test_promotion(a, p)
    end
    for p in (2x, 2x + 3)
        test_equal(LinearAlgebra.dot(a, p), adjoint(a) * p)
        test_equal(LinearAlgebra.dot(p, a), adjoint(p) * a)
    end
    p = MA.@rewrite((2x + 3) + a * (x + 1))
    @test MA.isequal_canonical(MP.coefficient(p, x), 2 + a)
    @test MA.isequal_canonical(MP.coefficient(p, one(x)), 3 + a)
    p = MA.@rewrite((2x + 3) - (x + 1) * a)
    @test MA.isequal_canonical(MP.coefficient(p, x), 2 - a)
    @test MA.isequal_canonical(MP.coefficient(p, one(x)), 3 - a)
    return
end

function test_scalar_equality(a, x)
    constant = a * one(x)
    for p in (constant, MP.polynomial(constant), (a * x) / x)
        for comp in (==, isequal)
            @test (@inferred comp(p, a))
            @test (@inferred comp(a, p))
        end
    end
    for p in (x, x^2, a * x, MP.polynomial(a * x), x / (x + 1))
        for comp in (==, isequal)
            @test !(@inferred comp(p, a))
            @test !(@inferred comp(a, p))
        end
    end
    # A variable index needs a coefficient type to construct its scalar zero.
    z = zero(1a)
    for p in (0 * x, zero(x), zero(x) / x), comp in (==, isequal)
        @test (@inferred comp(p, z))
        @test (@inferred comp(z, p))
    end
    return
end

allocated_product(a, b) = @allocated a * b

@testset "Registered scalar arithmetic" begin
    @test Symbol("@scalar_arithmetic") ∉ names(MP)
    @test isdefined(MP, Symbol("@scalar_arithmetic"))
    @test Base.get_extension(MP, :MultivariatePolynomialsJuMPExt) !== nothing
    @test Base.get_extension(
        MP,
        :MultivariatePolynomialsMathOptInterfaceExt,
    ) !== nothing
    model = JuMP.Model()
    JuMP.@variable(model, a)
    b = MOI.VariableIndex(1)
    DP.@polyvar x
    TP.@polyvar y
    for var in (x, y), scalar in (2, 1.5, 2 + im, big(2), a, 2a + 1, b, 2b + 1)
        @testset "$(typeof(var)), $(typeof(scalar))" begin
            test_scalar(scalar, var)
            test_scalar_equality(scalar, var)
            # MOI needs a coefficient type to add two variable indices.
            # Multiplication by 1 supplies it for the dot product's reduction.
            dot_scalar = scalar isa MOI.VariableIndex ? 1scalar : scalar
            test_dots(dot_scalar, [var, var])
            test_dots(dot_scalar, [var^2, var^3])
        end
    end
    for var in (x, y)
        @test_throws MethodError "unregistered scalar" * var
        @test_throws MethodError var + "unregistered scalar"
        for p in (var, var^2, 2var, var + 1, zero(var), var / (var + 1))
            @test promote_type(typeof(p), String) === Any
            @test p != UnregisteredMutable()
            @test UnregisteredMutable() != p
            @test !isequal(p, UnregisteredMutable())
            @test !isequal(UnregisteredMutable(), p)
        end
    end
    if VERSION >= v"1.13"
        for scalar in (2, 1.5), m in (y, y^2)
            allocated_product(scalar, m)
            allocated_product(m, scalar)
            @test allocated_product(scalar, m) == 0
            @test allocated_product(m, scalar) == 0
        end
    end
end

end

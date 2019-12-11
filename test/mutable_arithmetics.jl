import MutableArithmetics
const MA = MutableArithmetics

function all_tests(a, b, c, d, e, f, g)
    a_copy = deepcopy(a)
    b_copy = deepcopy(b)
    c_copy = deepcopy(c)
    d_copy = deepcopy(d)
    e_copy = deepcopy(e)
    f_copy = deepcopy(f)
    g_copy = deepcopy(g)
    @testset "Scalar" begin
        MA.Test.scalar_test(a)
        MA.Test.scalar_test(b)
        MA.Test.scalar_test(c)
        MA.Test.scalar_test(d)
        MA.Test.scalar_test(e)
        MA.Test.scalar_test(f)
        MA.Test.scalar_test(g)
    end
    @test isequal(a, a_copy)
    @test isequal(b, b_copy)
    @test isequal(c, c_copy)
    @test isequal(d, d_copy)
    @test isequal(e, e_copy)
    @testset "Quadratic" begin
        MA.Test.quadratic_test(a, b, c, d)
        MA.Test.quadratic_test(b, c, d, e)
        MA.Test.quadratic_test(c, d, e, f)
        MA.Test.quadratic_test(d, e, f, g)
        MA.Test.quadratic_test(e, f, g, a)
        MA.Test.quadratic_test(f, g, a, b)
        MA.Test.quadratic_test(g, a, b, c)
    end
    @test isequal(a, a_copy)
    @test isequal(b, b_copy)
    @test isequal(c, c_copy)
    @test isequal(d, d_copy)
    @test isequal(e, e_copy)
    @testset "Sparse" begin
        MA.Test.sparse_test(a, b, [a b c; b c a; a b a])
    end
    @test isequal(a, a_copy)
    @test isequal(b, b_copy)
    @test isequal(c, c_copy)
    @test isequal(d, d_copy)
    @test isequal(e, e_copy)
    @testset "Vector" begin
        MA.Test.array_test([a, b, c])
        MA.Test.array_test([b, c, d])
        MA.Test.array_test([c, d, e])
        MA.Test.array_test([d, e, a])
        MA.Test.array_test([e, a, b])
    end
    @test isequal(a, a_copy)
    @test isequal(b, b_copy)
    @test isequal(c, c_copy)
    @test isequal(d, d_copy)
    @test isequal(e, e_copy)
    @testset "Matrix" begin
        MA.Test.array_test([a b; c d])
        MA.Test.array_test([c e; e d])
        MA.Test.array_test([a b c; b c a; a b a])
        MA.Test.array_test([d b c; d c e; e b a])
    end
    @test isequal(a, a_copy)
    @test isequal(b, b_copy)
    @test isequal(c, c_copy)
    @test isequal(d, d_copy)
    @test isequal(e, e_copy)
end

Mod.@polyvar x y

@testset "MutableArithmetics with variables" begin
    # Creating 7 different variables here gives a long compile time for TypedPolynomials
    all_tests(x, y, x, y, x, y, x)
end

@testset "MutableArithmetics with monomials" begin
    a = x^2
    b = y^2
    c = x
    d = x^2 * y
    e = x^3
    f = x^5
    g = x^4
    all_tests(a, b, c, d, e, f, g)
end

@testset "MutableArithmetics with terms in $T" for T in [Int, BigInt]
    if MA.mutability(T) isa MA.IsMutable && MA.mutability(typeof(x * y)) isa MA.IsMutable
        @testset "Int" begin
            MA.Test.int_test(termtype(x, T), exclude = ["int_add", "int_add_mul", "int_zero"])
        end
    end
    a = T(2) * x^2
    b = T(4) * y^2
    c = T(3) * x
    d = T(5) * x^2 * y
    e = T(5) * x^3
    f = T(1) * x^5
    g = T(2) * x^4
    all_tests(a, b, c, d, e, f, g)
end

@testset "MutableArithmetics with polynomials in $T" for T in [Int, BigInt]
    if MA.mutability(T) isa MA.IsMutable
        @testset "Int" begin
            MA.Test.int_test(polynomialtype(x, T))
        end
    end
    a = T(2) * x^2 + T(3) * x * y + T(4) * y
    b = T(4) * y^2 - T(1) * x * y + T(4) * x
    c = T(1) * x^2 + T(3) * x - T(4)
    d = T(5) * x^2 * y - T(3) * x + T(4) * y
    e = T(5) * x^3 - T(3) * x + T(4)
    f = x^5 + x^4 + x^2 + x + T(1)
    g = T(2)x^4 + x^3
    all_tests(a, b, c, d, e, f, g)
end

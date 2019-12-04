import MutableArithmetics
const MA = MutableArithmetics

@testset "MutableArithmetics with $T" for T in [BigInt]
    Mod.@polyvar x y
#    @testset "Int" begin
        MA.Test.int_test(termtype(x, T), exclude = ["int_add", "int_add_mul", "int_zero"])
        MA.Test.int_test(polynomialtype(x, T))
#    end
#    a = T(2) * x^2 + T(3) * x * y + T(4) * y
#    b = T(4) * y^2 - T(1) * x * y + T(4) * x
#    c = T(1) * x^2 + T(3) * x - T(4)
#    d = T(5) * x^2 * y - T(3) * x + T(4) * y
#    @testset "Scalar" begin
#        MA.Test.scalar_test(a)
#    end
#    @testset "Quadratic" begin
#        MA.Test.quadratic_test(a, b, c, d)
#    end
#    @testset "Sparse" begin
#        MA.Test.sparse_test(a, b, T[a b c; b c a; a b a])
#    end
#    @testset "Vector" begin
#        MA.Test.array_test(T[a, b, c])
#    end
#    @testset "Matrix" begin
#        MA.Test.array_test(T[a b; c d])
#        MA.Test.array_test(T[a b c; b c a; a b a])
#    end
end

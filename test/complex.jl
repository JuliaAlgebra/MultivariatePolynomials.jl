@testset "ComplexPoly" begin
    Mod.@polyvar a
    @test isreal(a)
    @test !isrealpart(a)
    @test !isimagpart(a)
    @test !isconj(a)
    @test ordinary_variable(a) == a
    @test conj(a) == real(a) == a
    @test_broken iszero(imag(a)) # no iszero for MA.Zero
    @test isreal(a^3 + 5a^2 + 4a)
    @test !isreal(a^3 + 5im * a^2 + 4a)

    Mod.@complex_polyvar x y
    @test !isreal(x)
    @test !isrealpart(x) &&
          !isrealpart(conj(x)) &&
          isrealpart(real(x)) &&
          !isrealpart(imag(x))
    @test !isimagpart(x) &&
          !isimagpart(conj(x)) &&
          !isimagpart(real(x)) &&
          isimagpart(imag(x))
    @test !isconj(x) && isconj(conj(x)) && !isconj(real(x)) && !isconj(imag(x))
    @test ordinary_variable(conj(x)) ==
          ordinary_variable(real(x)) ==
          ordinary_variable(imag(x)) ==
          ordinary_variable(x)
    @testset "show" begin
        struct SymbolVar <: MP.AbstractVariable end
        struct SymbolConjVar <: MP.AbstractVariable end
        Base.isreal(::Union{SymbolVar,SymbolConjVar}) = false
        MP.name_base_indices(::Union{SymbolVar,SymbolConjVar}) = (:xy, (1, 2))
        MP.isconj(::SymbolConjVar) = true
        @test sprint(show, SymbolVar()) == "xy₁₋₂"
        @test sprint(show, SymbolConjVar()) == "x̅y̅₁₋₂"
    end
    @test conj(x) != x && conj(conj(x)) == x
    @test real(x) == real(conj(x))
    @test imag(conj(x)) == -imag(x)
    @test real(x * conj(x)) == real(x)^2 + imag(x)^2
    @test real(x + conj(x)) == 2real(x)
    @test iszero(imag(x + conj(x)))
    @test iszero(real(x - conj(x)))
    @test imag(x - conj(x)) == 2imag(x)
    @test !isreal(x^3 + 5x^2 + 4x)
    @test isreal(a^3 + 5a^2 + 4a + x^0)
    @test conj(7x^3 + 6y^2 - (3 + im)x + 8a) ==
          7conj(x)^3 + 6conj(y)^2 - (3 - im) * conj(x) + 8a
    @test conj(monomial_vector([x, y, x * y^3, y * x^4])) ==
          monomial_vector([conj(x), conj(y), conj(x * y^3), conj(y * x^4)])
    @test real(x * y^2 * a) ==
          (
        real(x) * real(y)^2 - imag(x) * 2 * real(y) * imag(y) -
        real(x) * imag(y)^2
    ) * a
    @test imag(x * y^2 * a) ==
          (
        imag(x) * real(y)^2 + real(x) * 2 * real(y) * imag(y) -
        imag(x) * imag(y)^2
    ) * a
    herm = [a x+y^2; conj(x)+conj(y)^2 x^2+conj(x)^2]
    @test herm == herm'

    @test degree_complex(x * y^2 * conj(y)^3) == 3
    @test degree_complex(x * y^2 * conj(y)^3, x) == 1
    @test degree_complex(x * y^2 * conj(y)^3, y) == 3
    @test degree_complex(a^5 * x * y^2 * conj(y)^4) == 9
    @test degree_complex(a^5 * x * y^2 * conj(y)^2) == 8
    @test halfdegree(x * y^2 * conj(y^3)) == 3
    @test halfdegree(x * a^5 * conj(y^2)) == 5
    @test halfdegree(x^2 * a^5 * conj(y^2)) == 5
    @test ordinary_variable([x, y, conj(x), a, real(x), imag(y)]) == [x, y, a]

    @test subs(4x + 8y^2 - 6x^3, [x, y] => [2 + 4im, 9 - im]) ==
          (1176 - 32im) * x^0
    @test_broken subs(4x + 8y^2 - 6x^3, [x, conj(y)] => [2 + 4im, 9 - im]) ==
          (1176 + 256im) * x^0
    @test_broken subs(4x + 8y^2 - 6x^3, [x, real(y)] => [2 + 4im, 9])
    @test_broken subs(4x + 8y^2 - 6x^3, [x, real(y)] => [2 + 4im, 9 + 0 * x^0]) ==
          1184 + 112im + 144im * imag(y) - 8imag(y)^2
end

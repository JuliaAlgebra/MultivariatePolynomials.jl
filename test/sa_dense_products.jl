module TestSADenseProducts

using Test
if !isdefined(Main, :Mod)
    include("SimplePolynomials.jl")
end
const Mod = isdefined(Main, :Mod) ? Main.Mod : SimplePolynomials
import MultivariatePolynomials as MP
import StarAlgebras as SA
import LinearAlgebra

@testset "Dense products use MA zero preparation" begin
    Mod.@polyvar x y z
    for T in (Int, BigInt)
        P = Mod.Polynomial{T}
        polys = P[x + 1 y + 2; z + 3 x + y]
        for (A, b) in (
            ([1 2; 3 4], P[x+1, y+2]),
            (polys, [2, 3]),
            (polys, P[y+1, z+2]),
            ([T(2)*x T(3)*y; T(4)*z T(1)*x], P[y+1, z+2]),
        )
            for B in (b, hcat(b, reverse(b), b))
                originals = MP.MA.mutable_copy.((A, B))
                expected = if B isa AbstractVector
                    [A[i, 1] * B[1] + A[i, 2] * B[2] for i in 1:2]
                else
                    [A[i, 1] * B[1, j] + A[i, 2] * B[2, j] for i in 1:2, j in 1:3]
                end
                result = @inferred A * B
                @test result == expected
                @test (@inferred MP.MA.operate(*, A, B)) == expected
                alg = parent(MP.MA.zero_for(eltype(result), A, B))
                @test all(r -> parent(r) == alg, result)
                output = similar(result)
                @test (@inferred LinearAlgebra.mul!(output, A, B)) === output
                @test output == expected
                scaled = similar(result)
                @test (@inferred LinearAlgebra.mul!(scaled, A, B, 2, 0)) ===
                      scaled
                @test scaled == 2 .* expected
                @test (@inferred MP.MA._mul!(scaled, A, B, 3)) === scaled
                @test scaled == 3 .* expected
                @test (A, B) == originals
                MP.MA.operate!(+, output[1], one(output[1]))
                @test output[2:end] == result[2:end]
                bad = B isa AbstractVector ? B[1:1] : B[1:1, :]
                original_output = MP.MA.mutable_copy(output)
                @test_throws DimensionMismatch LinearAlgebra.mul!(
                    output,
                    A,
                    bad,
                )
                @test output == original_output
            end
        end
        for B in (P[], Matrix{P}(undef, 0, 2))
            result = @inferred zeros(Int, 2, 0) * B
            @test all(iszero, result)
            MP.MA.operate!(+, result[1], one(result[1]))
            @test result[1] == 1
            @test all(iszero, result[2:end])
        end
        @test isempty(@inferred zeros(Int, 0, 2) * P[x+1, y+2])
        @test isempty(@inferred zeros(Int, 0, 2) * Matrix{P}(undef, 2, 0))
    end

    Mod.@complex_polyvar w
    A = Mod.Polynomial{Complex{Int}}[(1 + im)*w + x w + y; w + z (2 - im)*w + x]
    b = [1 + im, 2 - im]
    for left in (A, transpose(A), adjoint(A), view(A, :, :)),
        right in (b, hcat(b, reverse(b)))

        expected = if right isa AbstractVector
            [left[i, 1] * right[1] + left[i, 2] * right[2] for i in 1:2]
        else
            [
                left[i, 1] * right[1, j] + left[i, 2] * right[2, j] for
                i in 1:2, j in 1:2
            ]
        end
        @test (@inferred left * right) == expected
        output = similar(expected)
        @test (@inferred LinearAlgebra.mul!(output, left, right)) === output
        @test output == expected
    end
end

@testset "Outer products use MA zero preparation" begin
    Mod.@polyvar x y z
    for T in (Int, BigInt)
        P = Mod.Polynomial{T}
        left = P[x+1, y+2]
        for right in (transpose(P[y+1, z+2]), adjoint(P[y+1, z+2]))
            expected = left * right
            output = similar(expected)
            @test (@inferred MP.MA.operate_to!(output, *, left, right)) ===
                  output
            @test output == expected
            alg = parent(MP.MA.zero_for(P, left, right))
            @test all(p -> parent(p) == alg, output)
            MP.MA.operate!(+, output[1], one(output[1]))
            @test output[2:end] == expected[2:end]
        end
    end
end

end
